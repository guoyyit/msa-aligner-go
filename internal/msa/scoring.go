package msa

// SubstitutionMatrix is a 256x256 score table indexed by ASCII residue bytes.
type SubstitutionMatrix struct {
	Name   string
	Scores [256][256]float64
}

// Score returns the substitution score for residues a and b.
func (m *SubstitutionMatrix) Score(a, b byte) float64 {
	return m.Scores[a][b]
}

// aaOrder is the canonical BLOSUM amino-acid ordering.
var aaOrder = []byte{'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'}

// blosum62 is the standard BLOSUM62 matrix (half-bit units), symmetric.
var blosum62 = [20][20]float64{
	/* A */ {4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0},
	/* R */ {-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3},
	/* N */ {-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3},
	/* D */ {-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3},
	/* C */ {0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1},
	/* Q */ {-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2},
	/* E */ {-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2},
	/* G */ {0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3},
	/* H */ {-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3},
	/* I */ {-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3},
	/* L */ {-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1},
	/* K */ {-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2},
	/* M */ {-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1},
	/* F */ {-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1},
	/* P */ {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2},
	/* S */ {1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -2, -2, -1, 4, 1, -3, -2, -2},
	/* T */ {0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0},
	/* W */ {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3},
	/* Y */ {-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1},
	/* V */ {0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4},
}

// BLOSUM62 returns the standard BLOSUM62 substitution matrix.
func BLOSUM62() *SubstitutionMatrix {
	m := &SubstitutionMatrix{Name: "BLOSUM62"}
	idx := make(map[byte]int, len(aaOrder))
	for i, a := range aaOrder {
		idx[a] = i
	}
	for i, a := range aaOrder {
		for j, b := range aaOrder {
			m.Scores[a][b] = blosum62[i][j]
		}
	}
	return m
}

// dnaAmbiguity is the set of IUPAC nucleotide ambiguity codes.
var dnaAmbiguity = map[byte]bool{
	'R': true, 'Y': true, 'S': true, 'W': true, 'K': true, 'M': true,
	'B': true, 'D': true, 'H': true, 'V': true, 'N': true,
}

// DNAMatrix returns a nucleotide substitution matrix: equal residues score
// match, concrete mismatches score mismatch, and ambiguity/N pairs score 0.
func DNAMatrix(match, mismatch float64) *SubstitutionMatrix {
	m := &SubstitutionMatrix{Name: "DNA"}
	bases := []byte{'A', 'C', 'G', 'T', 'U'}
	for _, a := range bases {
		for _, b := range bases {
			if a == b {
				m.Scores[a][b] = match
			} else {
				m.Scores[a][b] = mismatch
			}
		}
	}
	// Ambiguity codes are treated as neutral (0) against concrete bases.
	for a := range dnaAmbiguity {
		for _, b := range bases {
			m.Scores[a][b] = 0
			m.Scores[b][a] = 0
		}
		m.Scores[a][a] = 0
	}
	m.Scores['N']['N'] = 0
	return m
}

// IdentityMatrix returns a match/mismatch matrix over all residues.
func IdentityMatrix(match, mismatch float64) *SubstitutionMatrix {
	m := &SubstitutionMatrix{Name: "identity"}
	for a := 0; a < 256; a++ {
		for b := 0; b < 256; b++ {
			if byte(a) == byte(b) {
				m.Scores[a][b] = match
			} else {
				m.Scores[a][b] = mismatch
			}
		}
	}
	return m
}
