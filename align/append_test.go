package align

import (
	"testing"

	"github.com/guoyyit/msa-aligner-go/cluster"
	"github.com/guoyyit/msa-aligner-go/distance"
	"github.com/guoyyit/msa-aligner-go/fasta"
)

func baseAlignment(t *testing.T) *Alignment {
	seqs := []fasta.Sequence{
		seq("a", "ACGTACGT"),
		seq("b", "ACGTTCGT"),
		seq("c", "ACGTACGA"),
	}
	m, err := distance.Compute(seqs, distance.Options{Concurrency: 1}, nil)
	if err != nil {
		t.Fatal(err)
	}
	tree, err := cluster.Build(m)
	if err != nil {
		t.Fatal(err)
	}
	a, err := Progressive(seqs, tree, DefaultParams(), nil)
	if err != nil {
		t.Fatal(err)
	}
	return a
}

func TestAppendPreservesExistingRowsByteForByte(t *testing.T) {
	base := baseAlignment(t)
	saved := map[string]string{}
	for i, id := range base.IDs {
		saved[id] = string(base.Rows[i])
	}
	out, err := Append(base, []fasta.Sequence{seq("new", "ACGTACGG")}, DefaultParams(), nil)
	if err != nil {
		t.Fatal(err)
	}
	for i, id := range out.IDs {
		if _, ok := saved[id]; !ok {
			continue
		}
		// Existing rows must be unchanged after stripping any appended
		// padding? No: they must not even gain padding gaps relative to
		// their residue order; but they may gain trailing gap columns.
		got := string(out.Rows[i])
		old := saved[id]
		if len(got) < len(old) {
			t.Fatalf("existing row %s shrank: %q -> %q", id, old, got)
		}
		// Prefix that predates the append must match byte-for-byte,
		// because Append never inserts columns into existing rows.
		if got[:len(old)] != old {
			t.Fatalf("existing row %s prefix changed: %q -> %q", id, old, got)
		}
		for k := len(old); k < len(got); k++ {
			if got[k] != '-' {
				t.Fatalf("existing row %s gained non-gap residue at new col %d", id, k)
			}
		}
	}
}

func TestAppendOnlyAddsRecordedColumns(t *testing.T) {
	base := baseAlignment(t)
	out, err := Append(base, []fasta.Sequence{seq("new", "TTTTACGT")}, DefaultParams(), nil)
	if err != nil {
		t.Fatal(err)
	}
	if out.Width < base.Width {
		t.Fatal("width must not shrink")
	}
	if len(out.AppendedColumns) != out.Width-base.Width {
		t.Fatalf("AppendedColumns=%v but width grew by %d", out.AppendedColumns, out.Width-base.Width)
	}
	for _, c := range out.AppendedColumns {
		if c < base.Width {
			t.Fatalf("appended column index %d points inside original width %d", c, base.Width)
		}
	}
}

func TestAppendMultipleAndWidths(t *testing.T) {
	base := baseAlignment(t)
	adds := []fasta.Sequence{
		seq("d", "GCGTACGT"),
		seq("e", "TTTTA"),
	}
	out, err := Append(base, adds, DefaultParams(), nil)
	if err != nil {
		t.Fatal(err)
	}
	if len(out.Rows) != len(base.Rows)+2 {
		t.Fatalf("row count = %d, want %d", len(out.Rows), len(base.Rows)+2)
	}
	for i, r := range out.Rows {
		if len(r) != out.Width {
			t.Fatalf("row %s width %d != %d", out.IDs[i], len(r), out.Width)
		}
		if ungap(r) != ungapExpected(out.IDs[i]) {
			t.Fatalf("residue order changed for %s: %q", out.IDs[i], r)
		}
	}
}

func ungapExpected(id string) string {
	data := map[string]string{
		"a": "ACGTACGT", "b": "ACGTTCGT", "c": "ACGTACGA",
		"d": "GCGTACGT", "e": "TTTTA",
	}
	return data[id]
}

func TestAppendDeterminismAndProgress(t *testing.T) {
	base := baseAlignment(t)
	adds := []fasta.Sequence{seq("d", "GCGTACGT"), seq("e", "TTTTACGT")}
	var done, total int
	out1, err := Append(base, adds, DefaultParams(), func(d, tot int) { done, total = d, tot })
	if err != nil {
		t.Fatal(err)
	}
	if done != 2 || total != 2 {
		t.Fatalf("progress = %d/%d", done, total)
	}
	out2, _ := Append(base, adds, DefaultParams(), nil)
	for i := range out1.Rows {
		if string(out1.Rows[i]) != string(out2.Rows[i]) {
			t.Fatalf("append not reproducible for row %s", out1.IDs[i])
		}
	}
}

func TestAppendDoesNotMutateInput(t *testing.T) {
	base := baseAlignment(t)
	before := map[string]string{}
	for i, id := range base.IDs {
		before[id] = string(base.Rows[i])
	}
	if _, err := Append(base, []fasta.Sequence{seq("new", "TTTTACGG")}, DefaultParams(), nil); err != nil {
		t.Fatal(err)
	}
	for i, id := range base.IDs {
		if string(base.Rows[i]) != before[id] {
			t.Fatalf("Append mutated existing alignment row %s", id)
		}
	}
}

func TestAppendValidation(t *testing.T) {
	base := baseAlignment(t)
	if _, err := Append(nil, []fasta.Sequence{seq("x", "ACGT")}, DefaultParams(), nil); err == nil {
		t.Fatal("nil existing must be rejected")
	}
	bad := Params{GapOpen: -1, GapExtend: 0}
	if _, err := Append(base, []fasta.Sequence{seq("x", "ACGT")}, bad, nil); err == nil {
		t.Fatal("missing substitution matrix must be rejected")
	}
}
