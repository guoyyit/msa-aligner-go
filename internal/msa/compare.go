package msa

import (
	"fmt"
	"sort"
	"strings"
)

// AlignmentDiff quantifies the bounded, explainable difference between an
// appended alignment and a full recomputation.
type AlignmentDiff struct {
	AppendScore    float64
	FullScore      float64
	ScoreDiff      float64
	ColumnsAdded   int // append length minus full length
	ColumnsShifted int // columns in append not found in full
	MaxShift       int // max column-index distance to a matching column
}

// Append aligns new sequences onto an existing alignment without modifying
// existing rows' residue-to-column relationships (gap columns may be inserted).
// New sequences are added one at a time via profile-single alignment.
func Append(existing *Alignment, newSeqs []Sequence, cfg Config, prog ProgressFunc) (*Alignment, error) {
	if cfg.Matrix == nil {
		cfg.Matrix = BLOSUM62()
	}
	if cfg.GapOpen == 0 {
		cfg.GapOpen = -10
	}
	if cfg.GapExtend == 0 {
		cfg.GapExtend = -1
	}
	if existing == nil || existing.Len == 0 {
		return nil, fmt.Errorf("msa: append requires a non-empty existing alignment")
	}

	rows := make([][]byte, len(existing.Seq))
	for i, s := range existing.Seq {
		rows[i] = []byte(s)
	}
	ids := append([]string{}, existing.IDs...)

	for k, ns := range newSeqs {
		merged, err := alignProfiles(rows, [][]byte{[]byte(ns.Seq)}, cfg)
		if err != nil {
			return nil, fmt.Errorf("msa: append %d: %w", k, err)
		}
		rows = merged
		ids = append(ids, ns.ID)
		if prog != nil {
			prog(Progress{
				Phase:   "append",
				Current: k + 1,
				Total:   len(newSeqs),
				Message: fmt.Sprintf("appended %d/%d", k+1, len(newSeqs)),
			})
		}
	}

	out := &Alignment{IDs: ids, Len: len(rows[0])}
	for _, r := range rows {
		if len(r) != out.Len {
			return nil, fmt.Errorf("msa: ragged append result")
		}
		out.Seq = append(out.Seq, string(r))
	}
	return out, nil
}

// CompareAlignments computes a bounded diff between two alignments that contain
// the same sequences. Rows are reordered to a canonical ID order before
// columns are compared.
func CompareAlignments(a, b *Alignment, cfg Config) (*AlignmentDiff, error) {
	aRows := canonicalRows(a)
	bRows := canonicalRows(b)
	if len(aRows) != len(bRows) {
		return nil, fmt.Errorf("msa: alignments have different sequence counts")
	}
	for i := range aRows {
		if aRows[i].id != bRows[i].id {
			return nil, fmt.Errorf("msa: alignments contain different sequences")
		}
	}

	aScore := AlignmentScore(a, cfg)
	bScore := AlignmentScore(b, cfg)

	aSigs := columnSignatures(aRows)
	bSigs := columnSignatures(bRows)

	bCount := make(map[string]int)
	for _, s := range bSigs {
		bCount[s]++
	}
	shifted := 0
	maxShift := 0
	for i, s := range aSigs {
		if bCount[s] > 0 {
			bCount[s]--
			continue
		}
		shifted++
		nearest := len(bSigs)
		for j, t := range bSigs {
			if t == s {
				if d := absInt(i - j); d < nearest {
					nearest = d
				}
			}
		}
		if nearest > maxShift {
			maxShift = nearest
		}
	}

	return &AlignmentDiff{
		AppendScore:    aScore,
		FullScore:      bScore,
		ScoreDiff:      aScore - bScore,
		ColumnsAdded:   a.Len - b.Len,
		ColumnsShifted: shifted,
		MaxShift:       maxShift,
	}, nil
}

// AlignmentScore computes the sum-of-pairs score of an alignment.
func AlignmentScore(a *Alignment, cfg Config) float64 {
	mat := cfg.Matrix
	if mat == nil {
		mat = BLOSUM62()
	}
	var total float64
	nRows := len(a.Seq)
	for c := 0; c < a.Len; c++ {
		for i := 0; i < nRows; i++ {
			for j := i + 1; j < nRows; j++ {
				x := a.Seq[i][c]
				y := a.Seq[j][c]
				if x == '-' || y == '-' {
					continue
				}
				total += mat.Score(x, y)
			}
		}
	}
	return total
}

type row struct {
	id  string
	seq string
}

func canonicalRows(a *Alignment) []row {
	rows := make([]row, len(a.IDs))
	for i := range a.IDs {
		rows[i] = row{id: a.IDs[i], seq: a.Seq[i]}
	}
	sort.Slice(rows, func(i, j int) bool { return rows[i].id < rows[j].id })
	return rows
}

func columnSignatures(rows []row) []string {
	if len(rows) == 0 {
		return nil
	}
	nCols := len(rows[0].seq)
	sigs := make([]string, nCols)
	for c := 0; c < nCols; c++ {
		var sb strings.Builder
		for _, r := range rows {
			sb.WriteByte(r.seq[c])
		}
		sigs[c] = sb.String()
	}
	return sigs
}

func absInt(x int) int {
	if x < 0 {
		return -x
	}
	return x
}
