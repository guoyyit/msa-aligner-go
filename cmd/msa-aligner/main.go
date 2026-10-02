// Command msa-aligner runs a reproducible multiple sequence alignment.
package main

import (
	"context"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"time"

	"github.com/guoyyit/msa-aligner-go/align"
	"github.com/guoyyit/msa-aligner-go/pipeline"
)

func main() {
	os.Exit(run(os.Args[1:]))
}

type flags struct {
	input       string
	batchDir    string
	alphabet    string
	k           int
	concurrency int
	match       int
	mismatch    int
	gapOpen     int
	gapExtend   int
	maxSeqs     int
	maxLen      int
	maxTotal    int64
	maxMem      int64
	checkpoint  string
	cpInterval  time.Duration
	output      string
	format      string
	quiet       bool
}

func parseFlags(args []string) (flags, *flag.FlagSet) {
	var f flags
	fs := flag.NewFlagSet("msa-aligner", flag.ContinueOnError)
	fs.StringVar(&f.input, "input", "-", "input FASTA file, or - for stdin")
	fs.StringVar(&f.batchDir, "batch", "", "directory of .fasta/.fa files to run as a batch")
	fs.StringVar(&f.alphabet, "alphabet", "dna", "alphabet: dna, rna, protein")
	fs.IntVar(&f.k, "k", 3, "k-mer length for pairwise distance")
	fs.IntVar(&f.concurrency, "concurrency", 1, "parallel distance workers")
	fs.IntVar(&f.match, "match", 2, "match score")
	fs.IntVar(&f.mismatch, "mismatch", -1, "mismatch score")
	fs.IntVar(&f.gapOpen, "gap-open", -10, "gap open penalty")
	fs.IntVar(&f.gapExtend, "gap-extend", -1, "gap extend penalty")
	fs.IntVar(&f.maxSeqs, "max-sequences", 5000, "maximum input sequence count")
	fs.IntVar(&f.maxLen, "max-length", 200000, "maximum single sequence length")
	fs.Int64Var(&f.maxTotal, "max-total-residues", 50_000_000, "maximum total residues")
	fs.Int64Var(&f.maxMem, "max-memory-bytes", 2*1024*1024*1024, "maximum distance matrix bytes")
	fs.StringVar(&f.checkpoint, "checkpoint", "", "checkpoint file enabling crash recovery")
	fs.DurationVar(&f.cpInterval, "checkpoint-interval", 5_000_000_000, "checkpoint write interval")
	fs.StringVar(&f.output, "output", "-", "output file, or - for stdout")
	fs.StringVar(&f.format, "format", "json", "output format: json or fasta")
	fs.BoolVar(&f.quiet, "quiet", false, "suppress progress on stderr")
	_ = fs.Parse(args) // parsed again by caller via the returned FlagSet
	return f, fs
}

func run(args []string) int {
	f, fs := parseFlags(args)
	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	cfg := pipeline.Config{
		Alphabet:    f.alphabet,
		K:           f.k,
		Concurrency: f.concurrency,
		Scoring: align.Params{
			Subst:     align.NewSubstMatrix(f.match, f.mismatch),
			GapOpen:   f.gapOpen,
			GapExtend: f.gapExtend,
		},
		Limits: pipeline.Limits{
			MaxSequences:      f.maxSeqs,
			MaxSequenceLength: f.maxLen,
			MaxTotalResidues:  f.maxTotal,
			MaxMemoryBytes:    f.maxMem,
		},
		CheckpointPath:     f.checkpoint,
		CheckpointInterval: f.cpInterval,
	}
	svc := pipeline.New(cfg)
	progress := progressSink(f.quiet)
	ctx := context.Background()

	if f.batchDir != "" {
		return runBatch(ctx, svc, f)
	}

	in, closeIn, err := openInput(f.input)
	if err != nil {
		fmt.Fprintln(os.Stderr, "error:", err)
		return 2
	}
	defer closeIn()
	res, err := svc.Align(ctx, in, progress)
	if err != nil {
		fmt.Fprintln(os.Stderr, "error:", err)
		return 1
	}
	if err := writeResult(f.output, f.format, res); err != nil {
		fmt.Fprintln(os.Stderr, "error:", err)
		return 1
	}
	return 0
}

func runBatch(ctx context.Context, svc *pipeline.Service, f flags) int {
	entries, err := os.ReadDir(f.batchDir)
	if err != nil {
		fmt.Fprintln(os.Stderr, "error:", err)
		return 2
	}
	var names []string
	for _, e := range entries {
		if e.IsDir() {
			continue
		}
		ext := strings.ToLower(filepath.Ext(e.Name()))
		if ext == ".fasta" || ext == ".fa" || ext == ".fna" {
			names = append(names, e.Name())
		}
	}
	sort.Strings(names)
	readers := map[string]io.Reader{}
	closers := map[string]*os.File{}
	for _, name := range names {
		fh, err := os.Open(filepath.Join(f.batchDir, name))
		if err != nil {
			fmt.Fprintln(os.Stderr, "error:", err)
			return 2
		}
		readers[name] = fh
		closers[name] = fh
	}
	defer func() {
		for _, fh := range closers {
			_ = fh.Close()
		}
	}()
	results, err := svc.AlignBatch(ctx, readers, progressSink(f.quiet))
	for _, name := range names {
		res, ok := results[name]
		if !ok {
			continue
		}
		base := strings.TrimSuffix(name, filepath.Ext(name))
		outName := base + ".aligned.json"
		if f.format == "fasta" {
			outName = base + ".aligned.fasta"
		}
		if err := writeResult(filepath.Join(f.batchDir, outName), f.format, res); err != nil {
			fmt.Fprintln(os.Stderr, "error writing", outName, ":", err)
			return 1
		}
	}
	if err != nil {
		fmt.Fprintln(os.Stderr, "batch completed with errors:", err)
		return 1
	}
	return 0
}

func openInput(path string) (io.Reader, func(), error) {
	if path == "-" {
		return os.Stdin, func() {}, nil
	}
	fh, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	return fh, func() { _ = fh.Close() }, nil
}

func progressSink(quiet bool) pipeline.ProgressFunc {
	if quiet {
		return nil
	}
	return func(p pipeline.Progress) {
		fmt.Fprintf(os.Stderr, "[%s] %d/%d %s\n", p.Stage, p.Done, p.Total, p.Detail)
	}
}

type jsonOutput struct {
	IDs          []string       `json:"ids"`
	Width        int            `json:"width"`
	Rows         []string       `json:"rows"`
	Consensus    string         `json:"consensus"`
	Columns      []columnOutput `json:"columns"`
	MergeOrder   []mergeOutput  `json:"merge_order"`
	AppendedCols []int          `json:"appended_columns,omitempty"`
}

type columnOutput struct {
	Index       int            `json:"index"`
	Consensus   string         `json:"consensus"`
	Identity    float64        `json:"identity"`
	Conserved   bool           `json:"conserved"`
	GapFraction float64        `json:"gap_fraction"`
	Counts      map[string]int `json:"counts"`
}

type mergeOutput struct {
	Step     int     `json:"step"`
	Left     int     `json:"left"`
	Right    int     `json:"right"`
	Distance float64 `json:"distance"`
}

func writeResult(path, format string, res *pipeline.Result) error {
	var data []byte
	var err error
	switch format {
	case "json":
		data, err = buildJSON(res)
	case "fasta":
		data = buildFASTA(res)
	default:
		return fmt.Errorf("unknown format %q", format)
	}
	if err != nil {
		return err
	}
	if path == "-" {
		_, err = os.Stdout.Write(data)
		return err
	}
	return os.WriteFile(path, data, 0o644)
}

func buildJSON(res *pipeline.Result) ([]byte, error) {
	a := res.Alignment
	out := jsonOutput{
		IDs:       a.IDs,
		Width:     a.Width,
		Rows:      make([]string, len(a.Rows)),
		Consensus: string(a.Consensus),
		Columns:   make([]columnOutput, len(a.Columns)),
	}
	for i, r := range a.Rows {
		out.Rows[i] = string(r)
	}
	for i, c := range a.Columns {
		counts := map[string]int{}
		for k, v := range c.Counts {
			counts[string(k)] = v
		}
		out.Columns[i] = columnOutput{
			Index: c.Index, Consensus: string(c.Consensus), Identity: c.Identity,
			Conserved: c.Conserved, GapFraction: c.GapFraction, Counts: counts,
		}
	}
	for _, s := range res.Tree.Steps {
		out.MergeOrder = append(out.MergeOrder, mergeOutput{
			Step: s.Step, Left: s.Left, Right: s.Right, Distance: s.Distance,
		})
	}
	out.AppendedCols = a.AppendedColumns
	return json.MarshalIndent(out, "", "  ")
}

func buildFASTA(res *pipeline.Result) []byte {
	var b strings.Builder
	a := res.Alignment
	for i, id := range a.IDs {
		fmt.Fprintf(&b, ">%s\n%s\n", id, a.Rows[i])
	}
	return []byte(b.String())
}
