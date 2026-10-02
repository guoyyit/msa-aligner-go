package msa

import "sort"

// ColumnConservation describes one alignment column.
type ColumnConservation struct {
	Position  int     // 1-based column index
	Consensus byte    // most common residue ('-' if column is all gaps)
	Conserved float64 // fraction of rows matching the consensus residue
	GapFrac   float64 // fraction of rows that are gaps
}

// Consensus returns the consensus sequence and per-column conservation.
// Ties for the most common residue are broken alphabetically for determinism.
func Consensus(a *Alignment) (string, []ColumnConservation) {
	if a == nil || a.Len == 0 {
		return "", nil
	}
	nRows := len(a.Seq)
	cols := make([]ColumnConservation, a.Len)
	consensus := make([]byte, a.Len)

	for c := 0; c < a.Len; c++ {
		counts := make(map[byte]int)
		gaps := 0
		for r := 0; r < nRows; r++ {
			ch := a.Seq[r][c]
			if ch == '-' {
				gaps++
			} else {
				counts[ch]++
			}
		}

		var best byte
		bestCount := -1
		residues := make([]byte, 0, len(counts))
		for res := range counts {
			residues = append(residues, res)
		}
		sort.Slice(residues, func(i, j int) bool { return residues[i] < residues[j] })
		for _, res := range residues {
			if counts[res] > bestCount {
				bestCount = counts[res]
				best = res
			}
		}

		cc := ColumnConservation{Position: c + 1, GapFrac: float64(gaps) / float64(nRows)}
		if bestCount <= 0 {
			cc.Consensus = '-'
			cc.Conserved = 0
		} else {
			cc.Consensus = best
			cc.Conserved = float64(bestCount) / float64(nRows)
		}
		cols[c] = cc
		consensus[c] = cc.Consensus
	}
	return string(consensus), cols
}
