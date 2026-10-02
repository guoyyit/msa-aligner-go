package pipeline

import (
	"context"
	"strings"
	"testing"
)

func TestPipelineAppendPreservesAndExplains(t *testing.T) {
	svc := New(DefaultConfig())
	res, err := svc.Align(context.Background(), strings.NewReader(fastaInput), nil)
	if err != nil {
		t.Fatal(err)
	}
	oldRows := map[string]string{}
	oldWidth := res.Alignment.Width
	for i, id := range res.Alignment.IDs {
		oldRows[id] = string(res.Alignment.Rows[i])
	}

	out, err := svc.AppendTo(context.Background(), res, strings.NewReader(">new1\nGCGTACGTACGT\n>new2\nTTTTAAAAGGGA\n"), nil)
	if err != nil {
		t.Fatal(err)
	}
	if len(out.Alignment.Rows) != len(res.Alignment.Rows)+2 {
		t.Fatalf("row count = %d", len(out.Alignment.Rows))
	}
	// Existing rows: prefix before oldWidth must be byte-identical.
	for _, id := range res.Alignment.IDs {
		var row string
		for i, rid := range out.Alignment.IDs {
			if rid == id {
				row = string(out.Alignment.Rows[i])
			}
		}
		if row[:oldWidth] != oldRows[id] {
			t.Fatalf("existing row %s changed within original width", id)
		}
	}
	// Every appended column index must point beyond the old width.
	for _, c := range out.Alignment.AppendedColumns {
		if c < oldWidth {
			t.Fatalf("appended column index %d inside old width %d", c, oldWidth)
		}
	}
	if out.Alignment.Width-oldWidth != len(out.Alignment.AppendedColumns) {
		t.Fatal("appended column count must explain the width difference")
	}
}

func TestPipelineAppendRejectsDuplicateID(t *testing.T) {
	svc := New(DefaultConfig())
	res, err := svc.Align(context.Background(), strings.NewReader(fastaInput), nil)
	if err != nil {
		t.Fatal(err)
	}
	// "a" already exists.
	if _, err := svc.AppendTo(context.Background(), res, strings.NewReader(">a\nGGGGGGGGGGGG\n"), nil); err == nil {
		t.Fatal("duplicate id in append must be rejected")
	}
	// Invalid character must also be rejected.
	if _, err := svc.AppendTo(context.Background(), res, strings.NewReader(">z\nZZZZZZZZZZZZ\n"), nil); err == nil {
		t.Fatal("illegal char in append must be rejected")
	}
}

func TestPipelineParameterChangeReproducible(t *testing.T) {
	makeCfg := func(gapOpen int) Config {
		c := DefaultConfig()
		c.Scoring.GapOpen = gapOpen
		return c
	}
	run := func(cfg Config) map[string]string {
		svc := New(cfg)
		res, err := svc.Align(context.Background(), strings.NewReader(fastaInput), nil)
		if err != nil {
			t.Fatal(err)
		}
		out := map[string]string{}
		for i, id := range res.Alignment.IDs {
			out[id] = string(res.Alignment.Rows[i])
		}
		return out
	}
	cfgA := makeCfg(-10)
	r1 := run(cfgA)
	r2 := run(cfgA)
	for id := range r1 {
		if r1[id] != r2[id] {
			t.Fatalf("same params produced different alignment for %s", id)
		}
	}

	cfgB := makeCfg(-2)
	r3 := run(cfgB)
	// Different gap penalties are allowed to change the result, but each
	// parameter set is individually reproducible (r3 vs a second r3).
	r4 := run(cfgB)
	for id := range r3 {
		if r3[id] != r4[id] {
			t.Fatalf("alternative params not reproducible for %s", id)
		}
	}
}
