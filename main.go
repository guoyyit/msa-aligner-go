package main

import (
	"context"
	"flag"
	"fmt"
	"io"
	"os"
	"runtime"

	"github.com/guoyyit/msa-aligner-go/internal/fasta"
	"github.com/guoyyit/msa-aligner-go/internal/msa"
)

func main() {
	if err := run(os.Args[1:]); err != nil {
		fmt.Fprintf(os.Stderr, "error: %v\n", err)
		os.Exit(1)
	}
}

func run(args []string) error {
	fs := flag.NewFlagSet("msa-aligner", flag.ContinueOnError)
	fs.SetOutput(os.Stderr)

	input := fs.String("i", "", "input FASTA file (required)")
	output := fs.String("o", "", "output aligned FASTA file (default: stdout)")
	consOut := fs.String("consensus", "", "consensus report output file (optional)")
	checkpointDir := fs.String("c", "", "checkpoint directory (enables crash recovery)")
	appendFile := fs.String("append", "", "existing alignment FASTA to append new sequences to")
	fullRecompute := fs.Bool("full", false, "also run full recompute and report diff (with -append)")

	matrixName := fs.String("matrix", "blosum62", "substitution matrix: blosum62 | dna | identity")
	gapOpen := fs.Float64("gap-open", -10, "gap open penalty")
	gapExtend := fs.Float64("gap-extend", -1, "gap extend penalty")
	kmer := fs.Int("k", 3, "k-mer size for pairwise distance")
	concurrency := fs.Int("concurrency", runtime.NumCPU(), "distance computation workers")

	maxSeqs := fs.Int("max-seqs", 10000, "maximum number of sequences")
	maxLen := fs.Int("max-len", 100000, "maximum sequence length")
	maxTotalRes := fs.Int64("max-total-res", 100_000_000, "maximum total residues")
	maxMemory := fs.Int64("max-memory", 1<<30, "estimated memory cap in bytes")

	if err := fs.Parse(args); err != nil {
		return err
	}
	if *input == "" {
		return fmt.Errorf("input file required (-i)")
	}

	data, err := os.ReadFile(*input)
	if err != nil {
		return fmt.Errorf("read input: %w", err)
	}
	seqs, err := fasta.Parse(data)
	if err != nil {
		return fmt.Errorf("parse input: %w", err)
	}

	cfg := msa.Config{
		GapOpen:      *gapOpen,
		GapExtend:    *gapExtend,
		KmerSize:     *kmer,
		Concurrency:  *concurrency,
		MaxSequences: *maxSeqs,
		MaxLength:    *maxLen,
		MaxTotalRes:  *maxTotalRes,
		MaxMemory:    *maxMemory,
	}
	switch *matrixName {
	case "blosum62":
		cfg.Matrix = msa.BLOSUM62()
	case "dna":
		cfg.Matrix = msa.DNAMatrix(2, -1)
	case "identity":
		cfg.Matrix = msa.IdentityMatrix(5, -4)
	default:
		return fmt.Errorf("unknown matrix %q (want blosum62|dna|identity)", *matrixName)
	}

	prog := func(p msa.Progress) {
		fmt.Fprintf(os.Stderr, "\r[%s] %s", p.Phase, p.Message)
	}

	var out *msa.Alignment
	var result *msa.Result

	if *appendFile != "" {
		existing, err := loadAlignment(*appendFile)
		if err != nil {
			return err
		}
		out, err = msa.Append(existing, seqs, cfg, prog)
		if err != nil {
			return err
		}
		if *fullRecompute {
			all := append(append([]msa.Sequence{}, toSequences(existing)...), seqs...)
			full, err := msa.Run(context.Background(), all, cfg, *checkpointDir, prog)
			if err != nil {
				return err
			}
			diff, err := msa.CompareAlignments(out, full.Alignment, cfg)
			if err != nil {
				return err
			}
			fmt.Fprintln(os.Stderr)
			fmt.Fprintf(os.Stderr, "diff vs full recompute: score %.2f vs %.2f (delta %.2f), shifted columns %d, max shift %d, length diff %d\n",
				diff.AppendScore, diff.FullScore, diff.ScoreDiff, diff.ColumnsShifted, diff.MaxShift, diff.ColumnsAdded)
		}
	} else {
		result, err = msa.Run(context.Background(), seqs, cfg, *checkpointDir, prog)
		if err != nil {
			return err
		}
		out = result.Alignment
	}

	fmt.Fprintln(os.Stderr)
	fmt.Fprintf(os.Stderr, "aligned %d sequences, %d columns\n", len(out.IDs), out.Len)

	w := os.Stdout
	if *output != "" {
		f, err := os.Create(*output)
		if err != nil {
			return fmt.Errorf("create output: %w", err)
		}
		defer f.Close()
		w = f
	}
	if err := msa.WriteFASTA(w, out); err != nil {
		return err
	}

	if *consOut != "" && result != nil {
		if err := writeConsensus(*consOut, result); err != nil {
			return err
		}
	}
	return nil
}

func loadAlignment(path string) (*msa.Alignment, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, fmt.Errorf("read alignment: %w", err)
	}
	seqs, err := fasta.ParseAlignment(data)
	if err != nil {
		return nil, fmt.Errorf("parse alignment: %w", err)
	}
	a := &msa.Alignment{IDs: make([]string, len(seqs)), Seq: make([]string, len(seqs))}
	for i, s := range seqs {
		a.IDs[i] = s.ID
		a.Seq[i] = s.Seq
	}
	if len(seqs) > 0 {
		a.Len = len(seqs[0].Seq)
	}
	return a, nil
}

func toSequences(a *msa.Alignment) []msa.Sequence {
	seqs := make([]msa.Sequence, len(a.IDs))
	for i := range a.IDs {
		// Strip gaps to recover raw residues for full recompute.
		raw := stripGaps(a.Seq[i])
		seqs[i] = msa.Sequence{ID: a.IDs[i], Seq: raw}
	}
	return seqs
}

func stripGaps(s string) string {
	var b []byte
	for i := 0; i < len(s); i++ {
		if s[i] != '-' && s[i] != '.' {
			b = append(b, s[i])
		}
	}
	return string(b)
}

func writeConsensus(path string, result *msa.Result) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	if _, err := fmt.Fprintf(f, "# consensus: %s\n", result.Consensus); err != nil {
		return err
	}
	if _, err := io.WriteString(f, "position\tconsensus\tconserved\tgap_fraction\n"); err != nil {
		return err
	}
	for _, c := range result.Columns {
		if _, err := fmt.Fprintf(f, "%d\t%c\t%.4f\t%.4f\n", c.Position, c.Consensus, c.Conserved, c.GapFrac); err != nil {
			return err
		}
	}
	return nil
}
