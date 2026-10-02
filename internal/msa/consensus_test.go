package msa

import "testing"

func TestConsensusFullyConserved(t *testing.T) {
	a := &Alignment{
		IDs: []string{"a", "b"},
		Seq: []string{"ACGT", "ACGT"},
		Len: 4,
	}
	cons, cols := Consensus(a)
	if cons != "ACGT" {
		t.Errorf("consensus = %q, want ACGT", cons)
	}
	for i, c := range cols {
		if c.Conserved != 1 {
			t.Errorf("col %d conserved = %v, want 1", i, c.Conserved)
		}
		if c.GapFrac != 0 {
			t.Errorf("col %d gapfrac = %v, want 0", i, c.GapFrac)
		}
	}
}

func TestConsensusMixed(t *testing.T) {
	a := &Alignment{
		IDs: []string{"a", "b", "c", "d"},
		Seq: []string{"A-", "AG", "A-", "AG"},
		Len: 2,
	}
	cons, cols := Consensus(a)
	if cols[0].Consensus != 'A' || cols[0].Conserved != 1 {
		t.Errorf("col0 = %+v", cols[0])
	}
	// col1: G,G,-,- -> consensus G, conserved 0.5, gapfrac 0.5
	if cols[1].Consensus != 'G' || cols[1].Conserved != 0.5 || cols[1].GapFrac != 0.5 {
		t.Errorf("col1 = %+v", cols[1])
	}
	if cons != "AG" {
		t.Errorf("consensus = %q, want AG", cons)
	}
}

func TestConsensusAllGap(t *testing.T) {
	a := &Alignment{
		IDs: []string{"a", "b"},
		Seq: []string{"-", "-"},
		Len: 1,
	}
	cons, cols := Consensus(a)
	if cols[0].Consensus != '-' || cols[0].Conserved != 0 {
		t.Errorf("col = %+v", cols[0])
	}
	if cons != "-" {
		t.Errorf("consensus = %q, want -", cons)
	}
}

func TestConsensusDeterministicTie(t *testing.T) {
	// Two residues tie; alphabetically first wins.
	a := &Alignment{
		IDs: []string{"a", "b"},
		Seq: []string{"C", "A"},
		Len: 1,
	}
	_, cols := Consensus(a)
	if cols[0].Consensus != 'A' {
		t.Errorf("consensus = %q, want A (tie-break)", cols[0].Consensus)
	}
}
