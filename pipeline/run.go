package pipeline

import (
	"context"
	"errors"
	"fmt"
	"io"
	"math"
	"os"
	"sync"
	"time"

	"github.com/guoyyit/msa-aligner-go/align"
	"github.com/guoyyit/msa-aligner-go/checkpoint"
	"github.com/guoyyit/msa-aligner-go/cluster"
	"github.com/guoyyit/msa-aligner-go/distance"
	"github.com/guoyyit/msa-aligner-go/fasta"
)

// ErrCheckpointMismatch means a checkpoint belongs to another job.
var ErrCheckpointMismatch = errors.New("pipeline: checkpoint input/config mismatch")

type runner struct {
	svc      *Service
	ctx      context.Context
	progress ProgressFunc
	jobID    string

	inputHash  string
	configHash string
	data       jobData
	restored   bool
	lastSave   time.Time
}

func (s *Service) run(ctx context.Context, input io.Reader, progress ProgressFunc) (*Result, error) {
	cfgHash, err := hashConfig(s.cfg)
	if err != nil {
		return nil, err
	}
	r := &runner{
		svc:        s,
		ctx:        ctx,
		progress:   progress,
		configHash: cfgHash,
	}

	r.report(stageParsing, 0, 1, "")
	parser := fasta.Parser{Alphabet: s.cfg.Alphabet}
	seqs, perr := parser.Parse(input)
	if perr != nil {
		return nil, perr
	}
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	// Parse always requires at least one valid record.
	if len(seqs) == 0 {
		return nil, errors.New("pipeline: no valid sequences")
	}
	r.inputHash = hashInput(seqs)
	r.jobID = "job-" + r.inputHash[:16]

	if err := checkLimits(seqs, s.cfg.Limits); err != nil {
		return nil, err
	}
	r.report(stageParsing, 1, 1, fmt.Sprintf("%d sequences accepted", len(seqs)))

	if s.cfg.CheckpointPath != "" {
		if err := r.restore(); err != nil {
			return nil, err
		}
	}

	matrix, err := r.computeMatrix(seqs)
	if err != nil {
		return nil, err
	}

	tree, skippedCluster := r.cachedTree(matrix)
	if !skippedCluster {
		r.report(stageCluster, 0, 1, "")
		tree, err = cluster.Build(matrix)
		if err != nil {
			return nil, err
		}
		r.data.Stage = stageCluster
		r.data.Steps = stepsToData(tree)
		r.data.Members = tree.Members
		if err := r.save(seqs, matrix, nil); err != nil {
			return nil, err
		}
		r.report(stageCluster, 1, 1, fmt.Sprintf("%d merges", len(tree.Steps)))
	}

	var alignment *align.Alignment
	if r.restored && r.data.Stage == stageDone {
		alignment = alignmentFromData(r.data)
		r.report(stageAlign, len(tree.Steps), len(tree.Steps), "resumed")
	} else {
		r.report(stageAlign, 0, max1(len(tree.Steps)), "")
		alignment, err = align.Progressive(seqs, tree, s.cfg.Scoring, func(done, total int) {
			r.report(stageAlign, done, total, "")
			r.maybeSave(seqs, matrix)
		})
		if err != nil {
			return nil, err
		}
		r.data.Stage = stageDone
		r.data.AlignedIDs = alignment.IDs
		r.data.Rows = rowsToStrings(alignment.Rows)
		r.data.Appended = alignment.AppendedColumns
		if err := r.save(seqs, matrix, alignment); err != nil {
			return nil, err
		}
		r.report(stageDone, 1, 1, "")
	}

	return &Result{
		IDs:       matrix.IDs,
		Sequences: seqs,
		Matrix:    matrix,
		Tree:      tree,
		Alignment: alignment,
	}, nil
}

func max1(x int) int {
	if x == 0 {
		return 1
	}
	return x
}

func (r *runner) report(stage string, done, total int, detail string) {
	if r.progress != nil {
		r.progress(Progress{Stage: stage, Done: done, Total: total, Detail: detail})
	}
}

func (r *runner) restore() error {
	st, err := checkpoint.Load(r.svc.cfg.CheckpointPath)
	if errors.Is(err, os.ErrNotExist) {
		return nil
	}
	if err != nil {
		return err // ErrCorrupt: refuse silent recovery
	}
	if st.InputHash != r.inputHash || st.ConfigHash != r.configHash {
		return fmt.Errorf("%w (stored input=%s, config=%s)", ErrCheckpointMismatch,
			short(st.InputHash), short(st.ConfigHash))
	}
	var d jobData
	if err := st.DecodeData(&d); err != nil {
		return err
	}
	r.data = d
	r.restored = true
	return nil
}

func short(s string) string {
	if len(s) > 8 {
		return s[:8]
	}
	return s
}

func (r *runner) cachedTree(m *distance.Matrix) (*cluster.Tree, bool) {
	if !r.restored || (r.data.Stage != stageCluster && r.data.Stage != stageAlign && r.data.Stage != stageDone) {
		return nil, false
	}
	if r.data.MatrixN != m.N || len(r.data.Steps) == 0 {
		return nil, false
	}
	return treeFromData(m.IDs, r.data), true
}

func (r *runner) computeMatrix(seqs []fasta.Sequence) (*distance.Matrix, error) {
	n := len(seqs)
	total := n * (n - 1) / 2

	if r.restored && r.data.MatrixN == n && len(r.data.Matrix) == n {
		m := matrixFromData(r.data)
		if allPairsKnown(m) {
			r.report(stageDistance, total, total, "resumed from checkpoint")
			return m, nil
		}
		r.report(stageDistance, countKnown(m), total, "resuming")
		return r.fillMissing(seqs, m, total)
	}

	// Fresh start.
	r.report(stageDistance, 0, total, "")
	r.data = newJobData(seqs)
	m := distance.NewMatrix(append([]string(nil), r.data.IDs...))
	return r.fillMissing(seqs, m, total)
}

func countKnown(m *distance.Matrix) int {
	n := m.N
	known := 0
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			_, failed := m.Failed(i, j)
			if failed || !math.IsNaN(m.At(i, j)) {
				known++
			}
		}
	}
	return known
}

func allPairsKnown(m *distance.Matrix) bool {
	return countKnown(m) == m.N*(m.N-1)/2
}

type pairResult struct {
	i, j int
	d    float64
	msg  string
}

func (r *runner) fillMissing(seqs []fasta.Sequence, m *distance.Matrix, total int) (*distance.Matrix, error) {
	n := len(seqs)
	type pair struct{ i, j int }
	var missing []pair
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			_, failed := m.Failed(i, j)
			if math.IsNaN(m.At(i, j)) && !failed {
				missing = append(missing, pair{i, j})
			}
		}
	}

	workers := r.svc.cfg.Concurrency
	if workers < 1 {
		workers = 1
	}
	tasks := make(chan pair)
	results := make(chan pairResult, len(missing))
	var wg sync.WaitGroup
	for w := 0; w < workers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for pr := range tasks {
				d, msg := distance.PairDistance(seqs[pr.i].Data, seqs[pr.j].Data, distance.Options{
					K: r.svc.cfg.K, Alphabet: r.svc.cfg.Alphabet,
				})
				select {
				case results <- pairResult{pr.i, pr.j, d, msg}:
				case <-r.ctx.Done():
					return
				}
			}
		}()
	}

	// Single dispatcher: all matrix writes happen on this goroutine, so
	// no mutex is needed between workers and checkpoint serialization.
	dispatchDone := make(chan struct{})
	go func() {
		defer close(dispatchDone)
		for _, pr := range missing {
			select {
			case tasks <- pr:
			case <-r.ctx.Done():
				return
			}
		}
	}()

	completed := countKnown(m)
	want := len(missing)
	received := 0
	var ticker *time.Ticker
	var tickCh <-chan time.Time
	if r.svc.cfg.CheckpointPath != "" && r.svc.cfg.CheckpointInterval > 0 {
		ticker = time.NewTicker(r.svc.cfg.CheckpointInterval)
		defer ticker.Stop()
		tickCh = ticker.C
	}

	var cancelled bool
	for received < want && !cancelled {
		select {
		case res := <-results:
			applyResult(m, r, res, n)
			completed++
			received++
			r.report(stageDistance, completed, total, "")
		case <-tickCh:
			_ = r.save(seqs, m, nil)
		case <-r.ctx.Done():
			cancelled = true
		}
	}
	<-dispatchDone
	close(tasks)
	wg.Wait()
	// Drain any results workers produced just before observing cancel.
drain:
	for {
		select {
		case res := <-results:
			applyResult(m, r, res, n)
			completed++
		default:
			break drain
		}
	}

	// Always persist completed pairs so resume never recomputes them.
	if err := r.save(seqs, m, nil); err != nil {
		return nil, err
	}
	if cancelled {
		return nil, r.ctx.Err()
	}
	return m, nil
}

func applyResult(m *distance.Matrix, r *runner, res pairResult, n int) {
	m.Cells[res.i*n+res.j] = res.d
	m.Cells[res.j*n+res.i] = res.d
	r.data.Matrix[res.i][res.j] = jsonFloat(res.d)
	r.data.Matrix[res.j][res.i] = jsonFloat(res.d)
	if res.msg != "" {
		m.Errors[[2]int{res.i, res.j}] = res.msg
		r.data.MatrixErrs[res.i][res.j] = res.msg
	}
}
