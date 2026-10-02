package msa

import (
	"context"
	"crypto/sha256"
	"encoding/hex"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"runtime"
)

// Result is the output of an alignment run.
type Result struct {
	Alignment  *Alignment
	Consensus  string
	Columns    []ColumnConservation
	Checkpoint string // path to latest checkpoint, if any
}

// Run performs the full MSA pipeline: validate inputs, compute distances,
// build merge order, align progressively, and report consensus. It checkpoints
// after each phase and each merge so it can resume without recomputation.
func Run(ctx context.Context, seqs []Sequence, cfg Config, checkpointDir string, prog ProgressFunc) (*Result, error) {
	if ctx == nil {
		ctx = context.Background()
	}
	if err := ValidateLimits(seqs, cfg); err != nil {
		return nil, err
	}
	if cfg.KmerSize == 0 {
		cfg.KmerSize = 3
	}
	if cfg.Concurrency == 0 {
		cfg.Concurrency = runtime.NumCPU()
	}
	if cfg.GapOpen == 0 {
		cfg.GapOpen = -10
	}
	if cfg.GapExtend == 0 {
		cfg.GapExtend = -1
	}
	if cfg.Matrix == nil {
		cfg.Matrix = BLOSUM62()
	}

	inputHash := InputHash(seqs)
	n := len(seqs)

	state := &Checkpoint{
		InputHash:   inputHash,
		KmerSize:    cfg.KmerSize,
		Concurrency: cfg.Concurrency,
		GapOpen:     cfg.GapOpen,
		GapExtend:   cfg.GapExtend,
		MatrixName:  cfg.Matrix.Name,
		SeqIDs:      make([]string, n),
		SeqData:     make([]string, n),
	}
	for i, s := range seqs {
		state.SeqIDs[i] = s.ID
		state.SeqData[i] = s.Seq
	}

	cpPath := ""
	if checkpointDir != "" {
		cpPath = filepath.Join(checkpointDir, "alignment.json")
		if _, err := os.Stat(cpPath); err == nil {
			loaded, lerr := LoadCheckpoint(cpPath)
			if lerr != nil {
				return nil, lerr // ErrCheckpointCorrupt
			}
			if loaded.InputHash != inputHash {
				return nil, fmt.Errorf("%w: input hash mismatch", ErrCheckpointInput)
			}
			if loaded.KmerSize == cfg.KmerSize &&
				loaded.GapOpen == cfg.GapOpen &&
				loaded.GapExtend == cfg.GapExtend &&
				loaded.MatrixName == cfg.Matrix.Name {
				state.DistanceDone = loaded.DistanceDone
				state.MergeDone = loaded.MergeDone
				state.NextMerge = loaded.NextMerge
				state.ClusterIDs = loaded.ClusterIDs
				state.ClusterSeq = loaded.ClusterSeq
				if state.DistanceDone {
					state.DistMatrix = loaded.DistMatrix
				}
				if state.MergeDone {
					state.Order = loaded.Order
				}
			}
			// If config differs, start fresh (recompute).
		}
	}

	if err := ctx.Err(); err != nil {
		return nil, err
	}

	// Phase 1: distances.
	if !state.DistanceDone {
		dm, errs := ComputeDistances(seqs, cfg.KmerSize, cfg.Concurrency, prog)
		if len(errs) > 0 {
			return nil, fmt.Errorf("msa: distance errors: %v", errs)
		}
		state.DistMatrix = dm.Dist
		state.DistanceDone = true
		if err := saveCheckpoint(cpPath, state); err != nil {
			return nil, err
		}
	}

	if err := ctx.Err(); err != nil {
		return nil, err
	}

	// Phase 2: merge order.
	if !state.MergeDone {
		dm := &DistanceMatrix{N: n, Dist: state.DistMatrix}
		order, err := BuildMergeOrder(dm)
		if err != nil {
			return nil, err
		}
		state.Order = order
		state.MergeDone = true
		if err := saveCheckpoint(cpPath, state); err != nil {
			return nil, err
		}
	}

	if err := ctx.Err(); err != nil {
		return nil, err
	}

	// Phase 3: progressive alignment (resumable per merge).
	if state.ClusterIDs == nil {
		state.ClusterIDs = make(map[int][]string)
		state.ClusterSeq = make(map[int][][]byte)
		for i := 0; i < n; i++ {
			state.ClusterIDs[i] = []string{seqs[i].ID}
			state.ClusterSeq[i] = [][]byte{[]byte(seqs[i].Seq)}
		}
	}
	for k := state.NextMerge; k < len(state.Order.Events); k++ {
		if err := ctx.Err(); err != nil {
			return nil, err
		}
		ev := state.Order.Events[k]
		aSeq := state.ClusterSeq[ev.Left]
		bSeq := state.ClusterSeq[ev.Right]
		if aSeq == nil || bSeq == nil {
			return nil, fmt.Errorf("msa: merge %d: missing cluster (%d or %d)", k, ev.Left, ev.Right)
		}
		merged, err := alignProfiles(aSeq, bSeq, cfg)
		if err != nil {
			return nil, fmt.Errorf("msa: merge %d: %w", k, err)
		}
		newID := n + k
		state.ClusterSeq[newID] = merged
		state.ClusterIDs[newID] = append(append([]string{}, state.ClusterIDs[ev.Left]...), state.ClusterIDs[ev.Right]...)
		state.NextMerge = k + 1
		if err := saveCheckpoint(cpPath, state); err != nil {
			return nil, err
		}
		if prog != nil {
			prog(Progress{
				Phase:   "merge",
				Current: k + 1,
				Total:   len(state.Order.Events),
				Message: fmt.Sprintf("merged %d/%d", k+1, len(state.Order.Events)),
			})
		}
	}

	finalID := 2*n - 2
	rows := state.ClusterSeq[finalID]
	ids := state.ClusterIDs[finalID]
	if rows == nil {
		return nil, fmt.Errorf("msa: final cluster %d missing", finalID)
	}
	out := &Alignment{IDs: ids, Len: len(rows[0])}
	for _, r := range rows {
		out.Seq = append(out.Seq, string(r))
	}

	cons, cols := Consensus(out)
	return &Result{
		Alignment:  out,
		Consensus:  cons,
		Columns:    cols,
		Checkpoint: cpPath,
	}, nil
}

func saveCheckpoint(path string, cp *Checkpoint) error {
	if path == "" {
		return nil
	}
	return SaveCheckpoint(path, cp)
}

// ValidateLimits checks resource limits before any work is done.
func ValidateLimits(seqs []Sequence, cfg Config) error {
	n := len(seqs)
	if cfg.MaxSequences > 0 && n > cfg.MaxSequences {
		return fmt.Errorf("%w: %d sequences > %d", ErrLimitExceeded, n, cfg.MaxSequences)
	}
	var totalRes int64
	maxLen := 0
	for _, s := range seqs {
		if len(s.Seq) > maxLen {
			maxLen = len(s.Seq)
		}
		totalRes += int64(len(s.Seq))
	}
	if cfg.MaxLength > 0 && maxLen > cfg.MaxLength {
		return fmt.Errorf("%w: max length %d > %d", ErrLimitExceeded, maxLen, cfg.MaxLength)
	}
	if cfg.MaxTotalRes > 0 && totalRes > cfg.MaxTotalRes {
		return fmt.Errorf("%w: total residues %d > %d", ErrLimitExceeded, totalRes, cfg.MaxTotalRes)
	}
	distBytes := int64(n) * int64(n) * 8
	profileBytes := (int64(n) * int64(n) / 2) * int64(maxLen)
	est := distBytes + profileBytes
	if cfg.MaxMemory > 0 && est > cfg.MaxMemory {
		return fmt.Errorf("%w: estimated memory %d bytes > %d", ErrLimitExceeded, est, cfg.MaxMemory)
	}
	return nil
}

// InputHash returns a stable SHA-256 hash of the input sequences.
func InputHash(seqs []Sequence) string {
	h := sha256.New()
	for _, s := range seqs {
		h.Write([]byte(s.ID))
		h.Write([]byte{0})
		h.Write([]byte(s.Seq))
		h.Write([]byte{0})
	}
	return hex.EncodeToString(h.Sum(nil))
}

// WriteFASTA writes an alignment in FASTA format with 60-column wrapping.
func WriteFASTA(w io.Writer, a *Alignment) error {
	for i, id := range a.IDs {
		if _, err := fmt.Fprintf(w, ">%s\n", id); err != nil {
			return err
		}
		seq := a.Seq[i]
		for len(seq) > 0 {
			chunk := seq
			if len(chunk) > 60 {
				chunk = chunk[:60]
			}
			if _, err := fmt.Fprintln(w, chunk); err != nil {
				return err
			}
			seq = seq[len(chunk):]
		}
	}
	return nil
}
