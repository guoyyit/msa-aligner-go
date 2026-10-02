package msa

import (
	"fmt"
	"math"
)

const (
	stateM = 0 // match: align A column with B column
	stateE = 1 // gap in B: A column inserted
	stateF = 2 // gap in A: B column inserted
)

// AlignProgressive aligns sequences following the given merge order using
// profile-profile alignment with affine gap penalties.
func AlignProgressive(seqs []Sequence, order *MergeOrder, cfg Config, prog ProgressFunc) (*Alignment, error) {
	n := len(seqs)
	if n == 0 {
		return nil, fmt.Errorf("msa: no sequences to align")
	}
	if cfg.Matrix == nil {
		cfg.Matrix = BLOSUM62()
	}
	if cfg.GapOpen == 0 {
		cfg.GapOpen = -10
	}
	if cfg.GapExtend == 0 {
		cfg.GapExtend = -1
	}

	clusterIDs := make(map[int][]string)
	clusterSeq := make(map[int][][]byte)
	for i := 0; i < n; i++ {
		clusterIDs[i] = []string{seqs[i].ID}
		clusterSeq[i] = [][]byte{[]byte(seqs[i].Seq)}
	}
	if n == 1 {
		return &Alignment{
			IDs: []string{seqs[0].ID},
			Seq: []string{seqs[0].Seq},
			Len: len(seqs[0].Seq),
		}, nil
	}

	for k, ev := range order.Events {
		aSeq := clusterSeq[ev.Left]
		bSeq := clusterSeq[ev.Right]
		if aSeq == nil || bSeq == nil {
			return nil, fmt.Errorf("msa: merge %d: missing cluster (%d or %d)", k, ev.Left, ev.Right)
		}
		merged, err := alignProfiles(aSeq, bSeq, cfg)
		if err != nil {
			return nil, fmt.Errorf("msa: merge %d: %w", k, err)
		}
		newID := n + k
		clusterSeq[newID] = merged
		aIDs := clusterIDs[ev.Left]
		bIDs := clusterIDs[ev.Right]
		clusterIDs[newID] = append(append([]string{}, aIDs...), bIDs...)
		if prog != nil {
			prog(Progress{
				Phase:   "merge",
				Current: k + 1,
				Total:   len(order.Events),
				Message: fmt.Sprintf("merged %d/%d", k+1, len(order.Events)),
			})
		}
	}

	finalID := 2*n - 2
	rows := clusterSeq[finalID]
	ids := clusterIDs[finalID]
	if rows == nil {
		return nil, fmt.Errorf("msa: final cluster %d missing", finalID)
	}
	out := &Alignment{IDs: ids, Len: len(rows[0])}
	for _, r := range rows {
		if len(r) != out.Len {
			return nil, fmt.Errorf("msa: ragged final alignment")
		}
		out.Seq = append(out.Seq, string(r))
	}
	return out, nil
}

// alignProfiles aligns two profiles (rows of equal length) via Gotoh DP with
// affine gap penalties. Profile-column scores average residue-residue
// substitution scores; gap-residue pairs in a matched column are neutral.
func alignProfiles(a, b [][]byte, cfg Config) ([][]byte, error) {
	nA := len(a)
	nB := len(b)
	if nA == 0 || nB == 0 {
		return nil, fmt.Errorf("msa: empty profile")
	}
	lenA := len(a[0])
	lenB := len(b[0])
	for _, row := range a {
		if len(row) != lenA {
			return nil, fmt.Errorf("msa: ragged profile")
		}
	}
	for _, row := range b {
		if len(row) != lenB {
			return nil, fmt.Errorf("msa: ragged profile")
		}
	}

	gapOpen := cfg.GapOpen
	gapExtend := cfg.GapExtend
	mat := cfg.Matrix
	if mat == nil {
		mat = BLOSUM62()
	}

	negInf := math.Inf(-1)
	M := make([][]float64, lenA+1)
	E := make([][]float64, lenA+1)
	F := make([][]float64, lenA+1)
	trM := make([][]byte, lenA+1)
	trE := make([][]byte, lenA+1)
	trF := make([][]byte, lenA+1)
	for i := 0; i <= lenA; i++ {
		M[i] = make([]float64, lenB+1)
		E[i] = make([]float64, lenB+1)
		F[i] = make([]float64, lenB+1)
		trM[i] = make([]byte, lenB+1)
		trE[i] = make([]byte, lenB+1)
		trF[i] = make([]byte, lenB+1)
	}

	colScore := func(i, j int) float64 {
		var sum float64
		var count int
		for r := 0; r < nA; r++ {
			ca := a[r][i]
			if ca == '-' {
				continue
			}
			for s := 0; s < nB; s++ {
				cb := b[s][j]
				if cb == '-' {
					continue
				}
				sum += mat.Score(ca, cb)
				count++
			}
		}
		if count == 0 {
			return 0
		}
		return sum / float64(count)
	}

	// Boundary conditions.
	M[0][0] = 0
	E[0][0] = negInf
	F[0][0] = negInf
	for i := 1; i <= lenA; i++ {
		M[i][0] = negInf
		F[i][0] = negInf
		E[i][0] = gapOpen + float64(i-1)*gapExtend
		trE[i][0] = stateE
	}
	for j := 1; j <= lenB; j++ {
		M[0][j] = negInf
		E[0][j] = negInf
		F[0][j] = gapOpen + float64(j-1)*gapExtend
		trF[0][j] = stateF
	}

	for i := 1; i <= lenA; i++ {
		for j := 1; j <= lenB; j++ {
			s := colScore(i-1, j-1)

			// Match: align A col i-1 with B col j-1.
			best := M[i-1][j-1]
			st := stateM
			if E[i-1][j-1] > best {
				best = E[i-1][j-1]
				st = stateE
			}
			if F[i-1][j-1] > best {
				best = F[i-1][j-1]
				st = stateF
			}
			M[i][j] = best + s
			trM[i][j] = byte(st)

			// Gap in B: A col i-1 inserted.
			bestE := E[i-1][j] + gapExtend
			stE := stateE
			if M[i-1][j]+gapOpen > bestE {
				bestE = M[i-1][j] + gapOpen
				stE = stateM
			}
			if F[i-1][j]+gapOpen > bestE {
				bestE = F[i-1][j] + gapOpen
				stE = stateF
			}
			E[i][j] = bestE
			trE[i][j] = byte(stE)

			// Gap in A: B col j-1 inserted.
			bestF := F[i][j-1] + gapExtend
			stF := stateF
			if M[i][j-1]+gapOpen > bestF {
				bestF = M[i][j-1] + gapOpen
				stF = stateM
			}
			if E[i][j-1]+gapOpen > bestF {
				bestF = E[i][j-1] + gapOpen
				stF = stateE
			}
			F[i][j] = bestF
			trF[i][j] = byte(stF)
		}
	}

	finalState := stateM
	finalScore := M[lenA][lenB]
	if E[lenA][lenB] > finalScore {
		finalScore = E[lenA][lenB]
		finalState = stateE
	}
	if F[lenA][lenB] > finalScore {
		finalScore = F[lenA][lenB]
		finalState = stateF
	}
	_ = finalScore

	// Traceback.
	var aCols, bCols []int
	i, j := lenA, lenB
	st := finalState
	for i > 0 || j > 0 {
		progressed := false
		switch st {
		case stateM:
			if i > 0 && j > 0 {
				aCols = append(aCols, i-1)
				bCols = append(bCols, j-1)
				st = int(trM[i][j])
				i--
				j--
				progressed = true
			}
		case stateE:
			if i > 0 {
				aCols = append(aCols, i-1)
				bCols = append(bCols, -1)
				st = int(trE[i][j])
				i--
				progressed = true
			}
		case stateF:
			if j > 0 {
				aCols = append(aCols, -1)
				bCols = append(bCols, j-1)
				st = int(trF[i][j])
				j--
				progressed = true
			}
		}
		if !progressed {
			break
		}
	}
	reverseInts(aCols)
	reverseInts(bCols)

	outA := make([][]byte, nA)
	outB := make([][]byte, nB)
	for r := 0; r < nA; r++ {
		row := make([]byte, 0, len(aCols))
		for _, c := range aCols {
			if c < 0 {
				row = append(row, '-')
			} else {
				row = append(row, a[r][c])
			}
		}
		outA[r] = row
	}
	for r := 0; r < nB; r++ {
		row := make([]byte, 0, len(bCols))
		for _, c := range bCols {
			if c < 0 {
				row = append(row, '-')
			} else {
				row = append(row, b[r][c])
			}
		}
		outB[r] = row
	}

	out := make([][]byte, 0, nA+nB)
	out = append(out, outA...)
	out = append(out, outB...)
	return out, nil
}

func reverseInts(s []int) {
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
}
