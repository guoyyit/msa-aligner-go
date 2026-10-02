package msa

import "testing"

func TestBLOSUM62(t *testing.T) {
	m := BLOSUM62()
	cases := []struct {
		a, b byte
		want float64
	}{
		{'A', 'A', 4},
		{'R', 'R', 5},
		{'W', 'W', 11},
		{'C', 'C', 9},
		{'A', 'R', -1},
		{'D', 'E', 2},
		{'H', 'Y', 2},
		{'I', 'L', 2},
		{'I', 'V', 3},
		{'F', 'Y', 3},
		{'K', 'R', 2},
		{'S', 'T', 1},
		{'M', 'L', 2},
	}
	for _, c := range cases {
		if got := m.Score(c.a, c.b); got != c.want {
			t.Errorf("Score(%c,%c) = %v, want %v", c.a, c.b, got, c.want)
		}
	}
	// Symmetric.
	if m.Score('A', 'R') != m.Score('R', 'A') {
		t.Error("matrix not symmetric")
	}
}

func TestDNAMatrix(t *testing.T) {
	m := DNAMatrix(2, -1)
	if m.Score('A', 'A') != 2 {
		t.Errorf("A/A = %v", m.Score('A', 'A'))
	}
	if m.Score('A', 'G') != -1 {
		t.Errorf("A/G = %v", m.Score('A', 'G'))
	}
	if m.Score('A', 'N') != 0 {
		t.Errorf("A/N = %v", m.Score('A', 'N'))
	}
}
