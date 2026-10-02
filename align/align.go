package align

import (
	"errors"

	"github.com/guoyyit/msa-aligner-go/cluster"
	"github.com/guoyyit/msa-aligner-go/fasta"
)

// SubstMatrix is a configurable residue-pair score table.
type SubstMatrix struct {
	scores   map[[2]byte]int
	match    int
	mismatch int
}

// NewSubstMatrix builds a table with a uniform match/mismatch score.
func NewSubstMatrix(match, mismatch int) *SubstMatrix {
	return &SubstMatrix{scores: map[[2]byte]int{}, match: match, mismatch: mismatch}
}

// Score returns the substitution score for a residue pair.
func (s *SubstMatrix) Score(a, b byte) int {
	if v, ok := s.scores[[2]byte{a, b}]; ok {
		return v
	}
	if v, ok := s.scores[[2]byte{b, a}]; ok {
		return v
	}
	if a == b {
		return s.match
	}
	return s.mismatch
}

// Set overrides the score for one residue pair.
func (s *SubstMatrix) Set(a, b byte, score int) {
	s.scores[[2]byte{a, b}] = score
}

// Equal reports whether two tables have identical behavior.
func (s *SubstMatrix) Equal(o *SubstMatrix) bool {
	if s == nil || o == nil {
		return s == o
	}
	if s.match != o.match || s.mismatch != o.mismatch || len(s.scores) != len(o.scores) {
		return false
	}
	for k, v := range s.scores {
		if ov, ok := o.scores[k]; !ok || ov != v {
			return false
		}
	}
	return true
}

// Params are the configurable alignment parameters.
type Params struct {
	Subst     *SubstMatrix
	GapOpen   int
	GapExtend int
}

// DefaultParams returns match +2 / mismatch -1, gap open -10 / extend -1.
func DefaultParams() Params {
	return Params{
		Subst:     NewSubstMatrix(2, -1),
		GapOpen:   -10,
		GapExtend: -1,
	}
}

// ColumnStat is the per-column consensus analysis.
type ColumnStat struct {
	Index       int
	Consensus   byte
	Identity    float64
	Conserved   bool
	GapFraction float64
	Counts      map[byte]int
}

// Alignment is a column-aligned set of sequences with column statistics.
type Alignment struct {
	IDs             []string
	Rows            [][]byte
	Width           int
	Consensus       []byte
	Columns         []ColumnStat
	AppendedColumns []int
}

// ProgressFunc reports completed/total alignment merges.
type ProgressFunc func(done, total int)

// ErrParams is returned when alignment parameters are unusable.
var ErrParams = errors.New("align: invalid parameters")

const negInf = -(1<<30 - 1)

// Pairwise aligns two raw sequences with affine-gap Needleman-Wunsch.
//
// States: diag (substitution), up (consume a, gap in b),
// left (consume b, gap in a). Traceback ties resolve in the fixed order
// diag > up > left, so identical input yields identical output everywhere.
func Pairwise(a, b []byte, p Params) (alignedA, alignedB []byte, score int) {
	n, m := len(a), len(b)
	diag, up, left := makeMatrix(n+1, m+1), makeMatrix(n+1, m+1), makeMatrix(n+1, m+1)
	diag[0][0], up[0][0], left[0][0] = 0, negInf, negInf
	for i := 1; i <= n; i++ {
		diag[i][0], up[i][0], left[i][0] = negInf, p.GapOpen+p.GapExtend*i, negInf
	}
	for j := 1; j <= m; j++ {
		diag[0][j], up[0][j], left[0][j] = negInf, negInf, p.GapOpen+p.GapExtend*j
	}
	for i := 1; i <= n; i++ {
		for j := 1; j <= m; j++ {
			sub := p.Subst.Score(a[i-1], b[j-1])
			diag[i][j] = max3(diag[i-1][j-1], up[i-1][j-1], left[i-1][j-1]) + sub
			up[i][j] = max2(diag[i-1][j]+p.GapOpen+p.GapExtend, up[i-1][j]+p.GapExtend)
			left[i][j] = max2(diag[i][j-1]+p.GapOpen+p.GapExtend, left[i][j-1]+p.GapExtend)
		}
	}
	end := max3(diag[n][m], up[n][m], left[n][m])

	// Predecessor state (0=diag, 1=up, 2=left) for every state cell.
	prevDiag := make([][]byte, n+1)
	prevUp := make([][]byte, n+1)
	prevLeft := make([][]byte, n+1)
	for i := range prevDiag {
		prevDiag[i] = make([]byte, m+1)
		prevUp[i] = make([]byte, m+1)
		prevLeft[i] = make([]byte, m+1)
	}
	for i := 1; i <= n; i++ {
		prevUp[i][0] = 1 // up(i,0) always extends the first gap run
	}
	for j := 1; j <= m; j++ {
		prevLeft[0][j] = 2
	}
	for i := 1; i <= n; i++ {
		for j := 1; j <= m; j++ {
			switch max3(diag[i-1][j-1], up[i-1][j-1], left[i-1][j-1]) {
			case diag[i-1][j-1]:
				prevDiag[i][j] = 0
			case up[i-1][j-1]:
				prevDiag[i][j] = 1
			default:
				prevDiag[i][j] = 2
			}
			if diag[i-1][j]+p.GapOpen+p.GapExtend >= up[i-1][j]+p.GapExtend {
				prevUp[i][j] = 0
			} else {
				prevUp[i][j] = 1
			}
			if diag[i][j-1]+p.GapOpen+p.GapExtend >= left[i][j-1]+p.GapExtend {
				prevLeft[i][j] = 0
			} else {
				prevLeft[i][j] = 2
			}
		}
	}

	i, j := n, m
	state := byte(0)
	switch end {
	case diag[n][m]:
		state = 0
	case up[n][m]:
		state = 1
	default:
		state = 2
	}
	var ra, rb []byte
	for i > 0 || j > 0 {
		if i > 0 && j > 0 {
			switch state {
			case 1:
				ra, rb = append(ra, a[i-1]), append(rb, '-')
				state = prevUp[i][j]
				i--
			case 2:
				ra, rb = append(ra, '-'), append(rb, b[j-1])
				state = prevLeft[i][j]
				j--
			default:
				ra, rb = append(ra, a[i-1]), append(rb, b[j-1])
				state = prevDiag[i][j]
				i--
				j--
			}
			continue
		}
		if i > 0 {
			ra, rb = append(ra, a[i-1]), append(rb, '-')
			i--
		} else {
			ra, rb = append(ra, '-'), append(rb, b[j-1])
			j--
		}
	}
	return reverse(ra), reverse(rb), end
}

func makeMatrix(n, m int) [][]int {
	rows := make([][]int, n)
	for i := range rows {
		rows[i] = make([]int, m)
	}
	return rows
}

func max2(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func max3(a, b, c int) int { return max2(max2(a, b), c) }

func reverse(b []byte) []byte {
	for l, r := 0, len(b)-1; l < r; l, r = l+1, r-1 {
		b[l], b[r] = b[r], b[l]
	}
	return b
}

// Progressive builds the full alignment following the merge tree.
func Progressive(seqs []fasta.Sequence, tree *cluster.Tree, p Params, progress ProgressFunc) (*Alignment, error) {
	if p.Subst == nil {
		return nil, ErrParams
	}
	n := tree.N
	if n == 0 || len(seqs) != n {
		return nil, errors.New("align: sequence count does not match tree leaves")
	}
	// nodeRows[node] holds the aligned rows for every subtree node.
	nodeRows := make([][][]byte, 2*n-1)
	nodeIDs := make([][]string, 2*n-1)
	for i := 0; i < n; i++ {
		nodeRows[i] = [][]byte{append([]byte(nil), seqs[i].Data...)}
		nodeIDs[i] = []string{seqs[i].ID}
	}
	for step, ms := range tree.Steps {
		mergedRows := alignProfiles(nodeRows[ms.Left], nodeRows[ms.Right], p)
		mergedIDs := make([]string, 0, len(nodeIDs[ms.Left])+len(nodeIDs[ms.Right]))
		mergedIDs = append(mergedIDs, nodeIDs[ms.Left]...)
		mergedIDs = append(mergedIDs, nodeIDs[ms.Right]...)
		nodeRows[ms.NewIndex] = mergedRows
		nodeIDs[ms.NewIndex] = mergedIDs
		if progress != nil {
			progress(step+1, len(tree.Steps))
		}
	}
	if n == 1 {
		return finalizeRows([]string{seqs[0].ID}, [][]byte{append([]byte(nil), seqs[0].Data...)}), nil
	}
	root := tree.Steps[len(tree.Steps)-1].NewIndex
	return finalizeRows(nodeIDs[root], nodeRows[root]), nil
}

// alignProfiles aligns two already-aligned row blocks column by column,
// applying affine gaps with fixed diag>up>left tie-breaking.
func alignProfiles(x, y [][]byte, p Params) [][]byte {
	nx, ny := len(x), len(y)
	wX, wY := len(x[0]), len(y[0])
	colScore := func(i, j int) int {
		s := 0
		for a := 0; a < nx; a++ {
			for b := 0; b < ny; b++ {
				ca, cb := x[a][i], y[b][j]
				if ca != '-' && cb != '-' {
					s += p.Subst.Score(ca, cb)
				}
			}
		}
		return s
	}
	n, m := wX, wY
	diag, up, left := makeMatrix(n+1, m+1), makeMatrix(n+1, m+1), makeMatrix(n+1, m+1)
	diag[0][0], up[0][0], left[0][0] = 0, negInf, negInf
	for i := 1; i <= n; i++ {
		diag[i][0], up[i][0], left[i][0] = negInf, (p.GapOpen+p.GapExtend*i)*ny, negInf
	}
	for j := 1; j <= m; j++ {
		diag[0][j], up[0][j], left[0][j] = negInf, negInf, (p.GapOpen+p.GapExtend*j)*nx
	}
	for i := 1; i <= n; i++ {
		for j := 1; j <= m; j++ {
			s := colScore(i-1, j-1)
			diag[i][j] = max3(diag[i-1][j-1], up[i-1][j-1], left[i-1][j-1]) + s
			up[i][j] = max2(diag[i-1][j]+(p.GapOpen+p.GapExtend)*ny, up[i-1][j]+p.GapExtend*ny)
			left[i][j] = max2(diag[i][j-1]+(p.GapOpen+p.GapExtend)*nx, left[i][j-1]+p.GapExtend*nx)
		}
	}

	prevDiag, prevUp, prevLeft := make([][]byte, n+1), make([][]byte, n+1), make([][]byte, n+1)
	for i := range prevDiag {
		prevDiag[i], prevUp[i], prevLeft[i] = make([]byte, m+1), make([]byte, m+1), make([]byte, m+1)
	}
	for i := 1; i <= n; i++ {
		prevUp[i][0] = 1
	}
	for j := 1; j <= m; j++ {
		prevLeft[0][j] = 2
	}
	for i := 1; i <= n; i++ {
		for j := 1; j <= m; j++ {
			switch max3(diag[i-1][j-1], up[i-1][j-1], left[i-1][j-1]) {
			case diag[i-1][j-1]:
				prevDiag[i][j] = 0
			case up[i-1][j-1]:
				prevDiag[i][j] = 1
			default:
				prevDiag[i][j] = 2
			}
			if diag[i-1][j]+(p.GapOpen+p.GapExtend)*ny >= up[i-1][j]+p.GapExtend*ny {
				prevUp[i][j] = 0
			} else {
				prevUp[i][j] = 1
			}
			if diag[i][j-1]+(p.GapOpen+p.GapExtend)*nx >= left[i][j-1]+p.GapExtend*nx {
				prevLeft[i][j] = 0
			} else {
				prevLeft[i][j] = 2
			}
		}
	}

	// Traceback yields two column masks describing how to pad each profile.
	var path []byte
	i, j := n, m
	end := max3(diag[n][m], up[n][m], left[n][m])
	state := byte(0)
	switch end {
	case up[n][m]:
		state = 1
	case left[n][m]:
		state = 2
	}
	for i > 0 || j > 0 {
		if i > 0 && j > 0 {
			path = append(path, state)
			switch state {
			case 1:
				i--
				state = prevUp[i+1][j]
			case 2:
				j--
				state = prevLeft[i][j+1]
			default:
				i--
				j--
				state = prevDiag[i+1][j+1]
			}
			continue
		}
		if i > 0 {
			path = append(path, 1)
			i--
		} else {
			path = append(path, 2)
			j--
		}
	}

	out := make([][]byte, 0, nx+ny)
	for r := 0; r < nx; r++ {
		out = append(out, make([]byte, 0, len(path)))
	}
	for r := 0; r < ny; r++ {
		out = append(out, make([]byte, 0, len(path)))
	}
	ci, cj := 0, 0
	for k := len(path) - 1; k >= 0; k-- {
		switch path[k] {
		case 0:
			for r := 0; r < nx; r++ {
				out[r] = append(out[r], x[r][ci])
			}
			for r := 0; r < ny; r++ {
				out[nx+r] = append(out[nx+r], y[r][cj])
			}
			ci++
			cj++
		case 1: // x column, y rows gap
			for r := 0; r < nx; r++ {
				out[r] = append(out[r], x[r][ci])
			}
			for r := 0; r < ny; r++ {
				out[nx+r] = append(out[nx+r], '-')
			}
			ci++
		case 2:
			for r := 0; r < nx; r++ {
				out[r] = append(out[r], '-')
			}
			for r := 0; r < ny; r++ {
				out[nx+r] = append(out[nx+r], y[r][cj])
			}
			cj++
		}
	}
	return out
}

func finalizeRows(ids []string, rows [][]byte) *Alignment {
	a := &Alignment{IDs: ids, Rows: rows}
	if len(rows) == 0 {
		return a
	}
	a.Width = len(rows[0])
	a.Columns = analyzeColumns(rows)
	a.Consensus = make([]byte, a.Width)
	for c := range a.Columns {
		a.Consensus[c] = a.Columns[c].Consensus
	}
	return a
}

// RestoreAlignment rebuilds an Alignment from previously aligned rows,
// recomputing width, consensus, and per-column statistics from the rows
// (which are the sole source of truth for derived data).
func RestoreAlignment(ids []string, rows [][]byte) *Alignment {
	return finalizeRows(append([]string(nil), ids...), rows)
}

// analyzeColumns produces the deterministic per-column consensus report.
// Consensus ties resolve by ascending residue byte; gap columns use '-'.
func analyzeColumns(rows [][]byte) []ColumnStat {
	if len(rows) == 0 {
		return nil
	}
	width := len(rows[0])
	cols := make([]ColumnStat, width)
	for c := 0; c < width; c++ {
		counts := map[byte]int{}
		gaps := 0
		for _, row := range rows {
			ch := row[c]
			counts[ch]++
			if ch == '-' {
				gaps++
			}
		}
		var consensus byte = 255
		best := -1
		for ch, n := range counts {
			if ch == '-' {
				continue
			}
			if n > best || (n == best && ch < consensus) {
				best, consensus = n, ch
			}
		}
		if best <= 0 {
			consensus = '-'
			best = gaps
		}
		identity := 0.0
		if consensus != '-' {
			identity = float64(best) / float64(len(rows))
		}
		cols[c] = ColumnStat{
			Index:       c,
			Consensus:   consensus,
			Identity:    identity,
			Conserved:   identity >= 1.0 && gaps == 0,
			GapFraction: float64(gaps) / float64(len(rows)),
			Counts:      counts,
		}
	}
	return cols
}

// Append merges additions (already in their final order) into an existing
// alignment without changing relative positions of existing residue columns.
//
// Semantics and the bounded-difference guarantee:
//   - Existing rows are copied byte-for-byte and never receive new gaps;
//     their residue order and column positions are exactly preserved.
//   - Each new sequence is aligned against the frozen consensus profile of
//     the existing alignment. Columns can only be added (never removed or
//     reordered), and Alignment.AppendedColumns records exactly which
//     columns are new, making the difference from a full recomputation
//     bounded and explainable.
//   - After the first addition the profile is frozen again, so additions
//     supplied in different orders each attach to the same fixed profile.
func Append(existing *Alignment, additions []fasta.Sequence, p Params, progress ProgressFunc) (*Alignment, error) {
	if existing == nil {
		return nil, errors.New("align: append requires an existing alignment")
	}
	if p.Subst == nil {
		return nil, ErrParams
	}
	if len(existing.Rows) == 0 || existing.Width == 0 {
		return nil, errors.New("align: cannot append to an empty alignment")
	}
	width := existing.Width
	consensus := append([]byte(nil), existing.Consensus...)
	if len(consensus) != width {
		consensus = make([]byte, width)
		for c := range consensus {
			consensus[c] = existing.Columns[c].Consensus
		}
	}

	// rows holds every output row; base rows first, then appended rows.
	rows := make([][]byte, 0, len(existing.Rows)+len(additions))
	for _, r := range existing.Rows {
		rows = append(rows, append([]byte(nil), r...))
	}
	ids := append([]string(nil), existing.IDs...)

	// mask[c] is true where output column c is an appended (new) column.
	mask := make([]bool, width)
	total := len(additions)
	for idx, add := range additions {
		aligned, alignedConsensus := alignToProfile(add.Data, consensus, p)
		newWidth := len(aligned)
		if newWidth != width {
			// Every existing row is aligned to the current consensus
			// (length == width), so each new ('-') marker simply adds a
			// gap; this holds for base rows and previously appended rows.
			for r := range rows {
				rows[r] = padAtNewColumns(rows[r], alignedConsensus)
			}
			mask = expandMask(mask, alignedConsensus)
			consensus = expandConsensus(consensus, alignedConsensus)
			width = newWidth
		}
		rows = append(rows, aligned)
		ids = append(ids, add.ID)
		if progress != nil {
			progress(idx+1, total)
		}
	}

	out := finalizeRows(ids, rows)
	for c, isNew := range mask {
		if isNew {
			out.AppendedColumns = append(out.AppendedColumns, c)
		}
	}
	return out, nil
}

// alignToProfile aligns seq against a fixed consensus profile using
// affine-gap DP. Returns the aligned sequence row and the aligned
// consensus row ('-' marks newly inserted columns).
func alignToProfile(seq, consensus []byte, p Params) (alignedSeq, alignedConsensus []byte) {
	n, m := len(seq), len(consensus)
	diag, up, left := makeMatrix(n+1, m+1), makeMatrix(n+1, m+1), makeMatrix(n+1, m+1)
	diag[0][0], up[0][0], left[0][0] = 0, negInf, negInf
	for i := 1; i <= n; i++ {
		diag[i][0], up[i][0], left[i][0] = negInf, p.GapOpen+p.GapExtend*i, negInf
	}
	for j := 1; j <= m; j++ {
		diag[0][j], up[0][j], left[0][j] = negInf, negInf, p.GapOpen+p.GapExtend*j
	}
	for i := 1; i <= n; i++ {
		for j := 1; j <= m; j++ {
			s := p.Subst.Score(seq[i-1], consensus[j-1])
			diag[i][j] = max3(diag[i-1][j-1], up[i-1][j-1], left[i-1][j-1]) + s
			up[i][j] = max2(diag[i-1][j]+p.GapOpen+p.GapExtend, up[i-1][j]+p.GapExtend)
			left[i][j] = max2(diag[i][j-1]+p.GapOpen+p.GapExtend, left[i][j-1]+p.GapExtend)
		}
	}

	prevDiag, prevUp, prevLeft := make([][]byte, n+1), make([][]byte, n+1), make([][]byte, n+1)
	for i := range prevDiag {
		prevDiag[i], prevUp[i], prevLeft[i] = make([]byte, m+1), make([]byte, m+1), make([]byte, m+1)
	}
	for i := 1; i <= n; i++ {
		prevUp[i][0] = 1
	}
	for j := 1; j <= m; j++ {
		prevLeft[0][j] = 2
	}
	for i := 1; i <= n; i++ {
		for j := 1; j <= m; j++ {
			switch max3(diag[i-1][j-1], up[i-1][j-1], left[i-1][j-1]) {
			case diag[i-1][j-1]:
				prevDiag[i][j] = 0
			case up[i-1][j-1]:
				prevDiag[i][j] = 1
			default:
				prevDiag[i][j] = 2
			}
			if diag[i-1][j]+p.GapOpen+p.GapExtend >= up[i-1][j]+p.GapExtend {
				prevUp[i][j] = 0
			} else {
				prevUp[i][j] = 1
			}
			if diag[i][j-1]+p.GapOpen+p.GapExtend >= left[i][j-1]+p.GapExtend {
				prevLeft[i][j] = 0
			} else {
				prevLeft[i][j] = 2
			}
		}
	}

	i, j := n, m
	state := byte(0)
	switch max3(diag[n][m], up[n][m], left[n][m]) {
	case up[n][m]:
		state = 1
	case left[n][m]:
		state = 2
	}
	var rs, rc []byte
	for i > 0 || j > 0 {
		if i > 0 && j > 0 {
			switch state {
			case 1:
				rs, rc = append(rs, seq[i-1]), append(rc, '-')
				i--
				state = prevUp[i+1][j]
			case 2:
				rs, rc = append(rs, '-'), append(rc, consensus[j-1])
				j--
				state = prevLeft[i][j+1]
			default:
				rs, rc = append(rs, seq[i-1]), append(rc, consensus[j-1])
				i--
				j--
				state = prevDiag[i+1][j+1]
			}
			continue
		}
		if i > 0 {
			rs, rc = append(rs, seq[i-1]), append(rc, '-')
			i--
		} else {
			rs, rc = append(rs, '-'), append(rc, consensus[j-1])
			j--
		}
	}
	return reverse(rs), reverse(rc)
}

// padAtNewColumns inserts a gap at every new column marker ('-') in the
// aligned consensus; existing characters keep their exact relative order.
func padAtNewColumns(row []byte, alignedConsensus []byte) []byte {
	out := make([]byte, 0, len(alignedConsensus))
	k := 0
	for _, c := range alignedConsensus {
		if c == '-' {
			out = append(out, '-')
			continue
		}
		out = append(out, row[k])
		k++
	}
	return out
}

func expandMask(mask []bool, alignedConsensus []byte) []bool {
	out := make([]bool, 0, len(alignedConsensus))
	k := 0
	for _, c := range alignedConsensus {
		if c == '-' {
			out = append(out, true)
			continue
		}
		out = append(out, mask[k])
		k++
	}
	return out
}

func expandConsensus(consensus, alignedConsensus []byte) []byte {
	out := make([]byte, 0, len(alignedConsensus))
	k := 0
	for _, c := range alignedConsensus {
		if c == '-' {
			out = append(out, '-')
			continue
		}
		out = append(out, consensus[k])
		k++
	}
	return out
}
