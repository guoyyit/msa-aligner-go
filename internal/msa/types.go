package msa

import (
	"errors"

	"github.com/guoyyit/msa-aligner-go/internal/fasta"
)

// Sequence is an input sequence with identifier and raw residues.
type Sequence = fasta.Sequence

// Alignment is a set of aligned sequence rows, all of equal length,
// containing gap characters ('-').
type Alignment struct {
	IDs []string
	Seq []string // aligned rows; len(Seq[i]) all equal
	Len int      // alignment length (number of columns)
}

// Progress reports long-running task progress to the caller.
type Progress struct {
	Phase   string // e.g. "distance", "merge", "checkpoint"
	Current int    // current unit index
	Total   int    // total units
	Message string // human-readable detail
}

// ProgressFunc receives progress updates during long-running operations.
type ProgressFunc func(Progress)

// Config controls alignment behavior and resource limits.
type Config struct {
	Matrix      *SubstitutionMatrix
	GapOpen     float64
	GapExtend   float64
	KmerSize    int
	Concurrency int

	// Resource limits; zero means use defaults.
	MaxSequences int
	MaxLength    int
	MaxTotalRes  int64
	MaxMemory    int64 // bytes
}

// DefaultConfig returns a configuration with sensible defaults.
func DefaultConfig() Config {
	return Config{}
}

// Distinguishable service-level errors.
var (
	ErrLimitExceeded      = errors.New("msa: resource limit exceeded")
	ErrCheckpointCorrupt  = errors.New("msa: checkpoint corrupted")
	ErrCheckpointInput    = errors.New("msa: checkpoint input mismatch")
	ErrIncompleteDistance = errors.New("msa: distance matrix incomplete")
)
