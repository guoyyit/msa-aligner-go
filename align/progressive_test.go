package align

import (
	"testing"

	"github.com/guoyyit/msa-aligner-go/cluster"
	"github.com/guoyyit/msa-aligner-go/distance"
	"github.com/guoyyit/msa-aligner-go/fasta"
)

func seq(id, data string) fasta.Sequence {
	return fasta.Sequence{ID: id, Data: []byte(data)}
}

func treeFor(t *testing.T, seqs []fasta.Sequence) *cluster.Tree {
	t.Helper()
	m, err := distance.Compute(seqs, distance.Options{Concurrency: 1}, nil)
	if err != nil {
		t.Fatal(err)
	}
	tree, err := cluster.Build(m)
	if err != nil {
		t.Fatal(err)
	}
	return tree
}

func validAlignment(t *testing.T, a *Alignment, seqs []fasta.Sequence) {
	t.Helper()
	if len(a.Rows) != len(seqs) || len(a.IDs) != len(seqs) {
		t.Fatalf("row/id count mismatch: %d/%d vs %d", len(a.Rows), len(a.IDs), len(seqs))
	}
	width := a.Width
	byID := map[string]string{}
	inputOrder := map[string]bool{}
	for _, s := range seqs {
		byID[s.ID] = string(s.Data)
		inputOrder[s.ID] = true
	}
	for i, id := range a.IDs {
		if !inputOrder[id] {
			t.Fatalf("unexpected id %s in alignment", id)
		}
		if len(a.Rows[i]) != width {
			t.Fatalf("row %s has width %d, want %d", id, len(a.Rows[i]), width)
		}
		if ungap(a.Rows[i]) != byID[id] {
			t.Fatalf("residues changed for %s: %q", id, a.Rows[i])
		}
	}
	if len(a.Consensus) != width || len(a.Columns) != width {
		t.Fatal("consensus/columns width mismatch")
	}
}

func TestProgressiveBasic(t *testing.T) {
	seqs := []fasta.Sequence{
		seq("a", "ACGTACGT"),
		seq("b", "ACGTTCGT"),
		seq("c", "TTTTAAAA"),
	}
	tree := treeFor(t, seqs)
	a, err := Progressive(seqs, tree, DefaultParams(), nil)
	if err != nil {
		t.Fatal(err)
	}
	validAlignment(t, a, seqs)
}

func TestProgressiveReproducible(t *testing.T) {
	seqs := []fasta.Sequence{
		seq("a", "GATTACA"),
		seq("b", "GCATGCT"),
		seq("c", "GATCACA"),
		seq("d", "TTAGCTA"),
	}
	tree := treeFor(t, seqs)
	first, err := Progressive(seqs, tree, DefaultParams(), nil)
	if err != nil {
		t.Fatal(err)
	}
	for run := 0; run < 3; run++ {
		got, err := Progressive(seqs, tree, DefaultParams(), nil)
		if err != nil {
			t.Fatal(err)
		}
		for i := range got.Rows {
			if string(got.Rows[i]) != string(first.Rows[i]) {
				t.Fatalf("row %d differs across runs", i)
			}
		}
	}
}

func TestProgressiveMergeOrderChangesAlignment(t *testing.T) {
	// Reordering the input records must not change the result: the tree
	// and alignment are keyed by distance, not slice position.
	seqs := []fasta.Sequence{
		seq("a", "ACGTACGT"),
		seq("b", "ACGTTCGT"),
		seq("c", "TTTTAAAA"),
	}
	tree := treeFor(t, seqs)
	a1, _ := Progressive(seqs, tree, DefaultParams(), nil)

	reversed := []fasta.Sequence{seqs[2], seqs[1], seqs[0]}
	tree2 := treeFor(t, reversed)
	a2, _ := Progressive(reversed, tree2, DefaultParams(), nil)

	got := map[string]string{}
	for i, id := range a2.IDs {
		got[id] = string(a2.Rows[i])
	}
	for i, id := range a1.IDs {
		if got[id] != string(a1.Rows[i]) {
			t.Fatalf("id %s alignment differs after reorder: %q vs %q", id, a1.Rows[i], got[id])
		}
	}
}

func TestProgressiveColumnStats(t *testing.T) {
	seqs := []fasta.Sequence{seq("a", "AAAAAAAA"), seq("b", "AAAAAAAA")}
	tree := treeFor(t, seqs)
	a, err := Progressive(seqs, tree, DefaultParams(), nil)
	if err != nil {
		t.Fatal(err)
	}
	for _, c := range a.Columns {
		if c.Consensus != 'A' || c.Identity != 1.0 || !c.Conserved || c.GapFraction != 0 {
			t.Fatalf("fully conserved column reported wrong: %+v", c)
		}
	}
}

func TestProgressiveSingleSequence(t *testing.T) {
	seqs := []fasta.Sequence{seq("only", "ACGT")}
	tree := treeFor(t, seqs)
	a, err := Progressive(seqs, tree, DefaultParams(), nil)
	if err != nil {
		t.Fatal(err)
	}
	if string(a.Rows[0]) != "ACGT" || a.Width != 4 {
		t.Fatalf("single sequence must pass through, got %+v", a.Rows)
	}
}

func TestProgressiveParamsValidation(t *testing.T) {
	seqs := []fasta.Sequence{seq("a", "AAAAAAAA"), seq("b", "AAAAAAAA")}
	tree := treeFor(t, seqs)
	if _, err := Progressive(seqs, tree, Params{GapOpen: -1, GapExtend: 0}, nil); err == nil {
		t.Fatal("missing substitution matrix must be rejected")
	}
}
