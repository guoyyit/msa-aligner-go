package pipeline

import (
	"encoding/json"
	"math"
	"time"

	"github.com/guoyyit/msa-aligner-go/align"
	"github.com/guoyyit/msa-aligner-go/checkpoint"
	"github.com/guoyyit/msa-aligner-go/cluster"
	"github.com/guoyyit/msa-aligner-go/distance"
	"github.com/guoyyit/msa-aligner-go/fasta"
)

func jsonMarshal(v any) (json.RawMessage, error) {
	return json.Marshal(v)
}

func newJobData(seqs []fasta.Sequence) jobData {
	n := len(seqs)
	d := jobData{
		Stage:      stageDistance,
		IDs:        make([]string, n),
		Data:       make([][]byte, n),
		Descs:      make([]string, n),
		MatrixN:    n,
		Matrix:     make([][]jsonFloat, n),
		MatrixErrs: make([][]string, n),
	}
	for i, s := range seqs {
		d.IDs[i] = s.ID
		d.Data[i] = append([]byte(nil), s.Data...)
		d.Descs[i] = s.Description
		d.Matrix[i] = make([]jsonFloat, n)
		d.MatrixErrs[i] = make([]string, n)
		for j := range d.Matrix[i] {
			d.Matrix[i][j] = jsonFloat(math.NaN())
		}
		d.Matrix[i][i] = 0
	}
	return d
}

func matrixFromData(d jobData) *distance.Matrix {
	n := d.MatrixN
	m := distance.NewMatrix(append([]string(nil), d.IDs...))
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			m.Cells[i*n+j] = float64(d.Matrix[i][j])
		}
		if i < len(d.MatrixErrs) {
			for j := i + 1; j < n; j++ {
				if msg := d.MatrixErrs[i][j]; msg != "" {
					m.Errors[[2]int{i, j}] = msg
				}
			}
		}
	}
	return m
}

func stepsToData(tree *cluster.Tree) []stepJSON {
	out := make([]stepJSON, len(tree.Steps))
	for i, s := range tree.Steps {
		out[i] = stepJSON{
			Step: s.Step, Left: s.Left, Right: s.Right,
			Distance: s.Distance, NewIndex: s.NewIndex,
		}
	}
	return out
}

func treeFromData(ids []string, d jobData) *cluster.Tree {
	steps := make([]cluster.MergeStep, len(d.Steps))
	for i, s := range d.Steps {
		steps[i] = cluster.MergeStep{
			Step: s.Step, Left: s.Left, Right: s.Right,
			Distance: s.Distance, NewIndex: s.NewIndex,
		}
	}
	members := make([][]int, len(d.Members))
	for i, mem := range d.Members {
		members[i] = append([]int(nil), mem...)
	}
	return &cluster.Tree{N: d.MatrixN, IDs: append([]string(nil), ids...), Steps: steps, Members: members}
}

func rowsToStrings(rows [][]byte) []string {
	out := make([]string, len(rows))
	for i, r := range rows {
		out[i] = string(r)
	}
	return out
}

func alignmentFromData(d jobData) *align.Alignment {
	rows := make([][]byte, len(d.Rows))
	for i, r := range d.Rows {
		rows[i] = []byte(r)
	}
	a := &align.Alignment{
		IDs:             append([]string(nil), d.AlignedIDs...),
		Rows:            rows,
		AppendedColumns: append([]int(nil), d.Appended...),
	}
	// Recompute consensus/column stats deterministically; they are derived
	// data and are not trusted from the checkpoint.
	rebuilt := align.RestoreAlignment(a.IDs, rows)
	rebuilt.AppendedColumns = a.AppendedColumns
	return rebuilt
}

// save writes the current state when checkpointing is configured.
// Callers must ensure exclusive access to m / r.data (the main goroutine
// owns it; the distance phase holds the worker mutex).
func (r *runner) save(seqs []fasta.Sequence, m *distance.Matrix, _ *align.Alignment) error {
	if r.svc.cfg.CheckpointPath == "" {
		return nil
	}
	n := len(seqs)
	if len(r.data.Matrix) != n {
		// Fresh run before matrix skeleton existed.
		r.data = newJobData(seqs)
	}
	if m != nil {
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				r.data.Matrix[i][j] = jsonFloat(m.Cells[i*n+j])
			}
		}
		for key, msg := range m.Errors {
			r.data.MatrixErrs[key[0]][key[1]] = msg
		}
	}
	st := &checkpoint.State{
		JobID:      r.jobID,
		InputHash:  r.inputHash,
		ConfigHash: r.configHash,
		UpdatedAt:  time.Now().UTC(),
	}
	raw, err := jsonMarshal(r.data)
	if err != nil {
		return err
	}
	st.Data = raw
	if err := checkpoint.Save(r.svc.cfg.CheckpointPath, st); err != nil {
		return err
	}
	r.lastSave = time.Now()
	return nil
}

// maybeSave persists completed work at most once per CheckpointInterval.
func (r *runner) maybeSave(seqs []fasta.Sequence, m *distance.Matrix) {
	if r.svc.cfg.CheckpointPath == "" {
		return
	}
	interval := r.svc.cfg.CheckpointInterval
	if interval <= 0 || time.Since(r.lastSave) < interval {
		return
	}
	_ = r.save(seqs, m, nil)
}
