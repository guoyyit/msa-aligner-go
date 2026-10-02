package pipeline

import (
	"context"
	"fmt"
	"io"
	"sort"

	"github.com/guoyyit/msa-aligner-go/align"
	"github.com/guoyyit/msa-aligner-go/fasta"
)

// AlignBatch runs each named input independently. Jobs are executed in
// ascending name order so batch results and progress are reproducible.
// A failing job does not abort the batch: its error is collected, and
// successful results are still returned.
func (s *Service) AlignBatch(ctx context.Context, inputs map[string]io.Reader, progress ProgressFunc) (map[string]*Result, error) {
	names := make([]string, 0, len(inputs))
	for name := range inputs {
		names = append(names, name)
	}
	sort.Strings(names)

	results := make(map[string]*Result, len(names))
	var batchErrs BatchErrors
	for idx, name := range names {
		if err := ctx.Err(); err != nil {
			return results, err
		}
		pf := func(p Progress) {
			if progress != nil {
				progress(Progress{
					Stage:  p.Stage,
					Done:   p.Done,
					Total:  p.Total,
					Detail: fmt.Sprintf("batch %d/%d job=%s %s", idx+1, len(names), name, p.Detail),
				})
			}
		}
		res, err := s.Align(ctx, inputs[name], pf)
		if err != nil {
			batchErrs = append(batchErrs, JobError{Job: name, Err: err.Error()})
			continue
		}
		results[name] = res
	}
	if len(batchErrs) > 0 {
		return results, batchErrs
	}
	return results, nil
}

// JobError is one failed batch job.
type JobError struct {
	Job string `json:"job"`
	Err string `json:"error"`
}

// BatchErrors collects per-job failures without aborting other jobs.
type BatchErrors []JobError

func (e BatchErrors) Error() string {
	out := fmt.Sprintf("%d of batch jobs failed:", len(e))
	for _, je := range e {
		out += "\n  - " + je.Job + ": " + je.Err
	}
	return out
}

// AppendTo extends a previous result with additional FASTA sequences.
//
// New sequences are parsed and validated against the same alphabet and
// limits. Their alignment attaches to the frozen existing alignment via
// align.Append, so existing rows never gain inserted columns and the
// difference is exactly AppendedColumns (bounded and explainable).
// New-sequence duplicate IDs, including duplicates of existing IDs, are
// rejected before any alignment work.
func (s *Service) AppendTo(ctx context.Context, prev *Result, newInput io.Reader, progress ProgressFunc) (*Result, error) {
	if prev == nil || prev.Alignment == nil {
		return nil, fmt.Errorf("pipeline: append requires a completed previous result")
	}
	parser := fasta.Parser{Alphabet: s.cfg.Alphabet}
	adds, err := parser.Parse(newInput)
	if err != nil {
		return nil, err
	}
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	// Reject duplicate IDs against the existing alignment and within adds.
	known := make(map[string]bool, len(prev.IDs))
	for _, id := range prev.IDs {
		known[id] = true
	}
	for _, a := range adds {
		if known[a.ID] {
			return nil, &appendIDError{id: a.ID}
		}
		known[a.ID] = true
	}
	if err := checkLimits(appendSeqs(prev.Sequences, adds), s.cfg.Limits); err != nil {
		return nil, err
	}

	total := len(adds)
	merged, err := align.Append(prev.Alignment, adds, s.cfg.Scoring, func(done, t int) {
		if progress != nil {
			progress(Progress{Stage: stageAlign, Done: done, Total: total, Detail: "appending"})
		}
	})
	if err != nil {
		return nil, err
	}

	allSeqs := appendSeqs(prev.Sequences, adds)
	return &Result{
		IDs:       merged.IDs,
		Sequences: allSeqs,
		Matrix:    prev.Matrix,
		Tree:      prev.Tree,
		Alignment: merged,
	}, nil
}

type appendIDError struct{ id string }

func (e *appendIDError) Error() string {
	return fmt.Sprintf("pipeline: appended sequence duplicates existing or repeated id %q", e.id)
}

func appendSeqs(existing []fasta.Sequence, adds []fasta.Sequence) []fasta.Sequence {
	out := make([]fasta.Sequence, 0, len(existing)+len(adds))
	out = append(out, existing...)
	out = append(out, adds...)
	return out
}
