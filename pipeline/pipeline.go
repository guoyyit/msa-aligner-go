// Package pipeline orchestrates the full MSA pipeline: FASTA parsing,
// bounded-concurrency distance matrix, deterministic clustering,
// progressive alignment, checkpointed resume, and batch execution.
package pipeline

import (
	"context"
	"io"
	"time"

	"github.com/guoyyit/msa-aligner-go/align"
	"github.com/guoyyit/msa-aligner-go/cluster"
	"github.com/guoyyit/msa-aligner-go/distance"
	"github.com/guoyyit/msa-aligner-go/fasta"
)

// Limits are pre-checked resource ceilings.
type Limits struct {
	MaxSequences      int
	MaxSequenceLength int
	MaxTotalResidues  int64
	MaxMemoryBytes    int64
}

// Config is the complete reproducible job identity for outputs.
type Config struct {
	Alphabet           string
	K                  int
	Concurrency        int
	Scoring            align.Params
	Warnings           []string
	Limits             Limits
	CheckpointPath     string
	CheckpointInterval time.Duration
}

// Progress is one observable progress event.
type Progress struct {
	Stage  string
	Done   int
	Total  int
	Detail string
}

// ProgressFunc receives progress events.
type ProgressFunc func(Progress)

// Result is the complete, reproducible alignment result.
type Result struct {
	IDs       []string
	Sequences []fasta.Sequence
	Matrix    *distance.Matrix
	Tree      *cluster.Tree
	Alignment *align.Alignment
	Warnings  []string
}

// Service runs alignment jobs.
type Service struct {
	cfg Config
}

// New constructs a Service with defaults applied.
func New(cfg Config) *Service {
	return &Service{cfg: cfg.withDefaults()}
}

// Align runs one FASTA input to completion, resuming from a valid
// checkpoint when one is present.
func (s *Service) Align(ctx context.Context, input io.Reader, progress ProgressFunc) (*Result, error) {
	return s.run(ctx, input, progress)
}
