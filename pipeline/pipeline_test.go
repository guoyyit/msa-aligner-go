package pipeline

import (
	"context"
	"io"
	"strings"
	"testing"
)

const fastaInput = ">a\nACGTACGTACGT\n>b\nACGTTCGTACGT\n>c\nTTTTAAAAGGGG\n>d\nACGTACGTACGA\n"

func runAlign(t *testing.T, cfg Config, input string) (*Result, []Progress) {
	t.Helper()
	var events []Progress
	svc := New(cfg)
	res, err := svc.Align(context.Background(), strings.NewReader(input), func(p Progress) {
		events = append(events, p)
	})
	if err != nil {
		t.Fatalf("align failed: %v", err)
	}
	return res, events
}

func TestPipelineBasic(t *testing.T) {
	res, events := runAlign(t, DefaultConfig(), fastaInput)
	if len(res.Sequences) != 4 {
		t.Fatalf("want 4 seqs, got %d", len(res.Sequences))
	}
	if len(res.Alignment.Rows) != 4 || res.Alignment.Width == 0 {
		t.Fatalf("bad alignment: %+v", res.Alignment)
	}
	stages := map[string]bool{}
	for _, e := range events {
		stages[e.Stage] = true
	}
	for _, want := range []string{stageParsing, stageDistance, stageCluster, stageAlign, stageDone} {
		if !stages[want] {
			t.Fatalf("missing progress stage %q", want)
		}
	}
}

func TestPipelineReproducibleAcrossConcurrency(t *testing.T) {
	var refRows map[string]string
	for _, workers := range []int{1, 4} {
		cfg := DefaultConfig()
		cfg.Concurrency = workers
		res, _ := runAlign(t, cfg, fastaInput)
		rows := map[string]string{}
		for i, id := range res.Alignment.IDs {
			rows[id] = string(res.Alignment.Rows[i])
		}
		if refRows == nil {
			refRows = rows
			continue
		}
		if len(rows) != len(refRows) {
			t.Fatal("row count differs")
		}
		for id, r := range rows {
			if refRows[id] != r {
				t.Fatalf("alignment for %s differs with concurrency %d: %q vs %q", id, workers, r, refRows[id])
			}
		}
	}
}

func TestPipelineMergeOrderDeterministic(t *testing.T) {
	res1, _ := runAlign(t, DefaultConfig(), fastaInput)
	res2, _ := runAlign(t, DefaultConfig(), fastaInput)
	if len(res1.Tree.Steps) != len(res2.Tree.Steps) {
		t.Fatal("merge step count differs")
	}
	for i := range res1.Tree.Steps {
		if res1.Tree.Steps[i] != res2.Tree.Steps[i] {
			t.Fatalf("merge step %d differs: %+v vs %+v", i, res1.Tree.Steps[i], res2.Tree.Steps[i])
		}
	}
}

func TestPipelineRejectsBadInput(t *testing.T) {
	svc := New(DefaultConfig())
	cases := map[string]string{
		"empty sequence": ">a\n>b\nACGT\n",
		"illegal char":   ">a\nACGZ\n",
		"duplicate id":   ">a\nACGT\n>a\nTTTT\n",
		"missing id":     ">\nACGT\n",
	}
	for name, in := range cases {
		if _, err := svc.Align(context.Background(), strings.NewReader(in), nil); err == nil {
			t.Fatalf("%s: expected error", name)
		}
	}
}

func TestPipelineRejectsBeforeLimit(t *testing.T) {
	cfg := DefaultConfig()
	cfg.Limits.MaxSequences = 2
	svc := New(cfg)
	_, err := svc.Align(context.Background(), strings.NewReader(fastaInput), nil)
	if err == nil {
		t.Fatal("expected limit error")
	}
	if le, ok := err.(*ErrLimitExceeded); !ok {
		t.Fatalf("want *ErrLimitExceeded, got %T: %v", err, err)
	} else if !strings.Contains(le.What, "count") {
		t.Fatalf("unexpected limit target: %s", le.What)
	}
}

func TestPipelineRejectsLongSequence(t *testing.T) {
	cfg := DefaultConfig()
	cfg.Limits.MaxSequenceLength = 4
	svc := New(cfg)
	_, err := svc.Align(context.Background(), strings.NewReader(">a\nACGTACGT\n>b\nTTTTAAAA\n"), nil)
	if err == nil {
		t.Fatal("expected per-sequence length limit error")
	}
}

func TestPipelineBatch(t *testing.T) {
	svc := New(DefaultConfig())
	inputs := map[string]string{
		"job-z": fastaInput,
		"job-a": ">x\nAAAAAAAAAAAA\n>y\nAAAATAAAAAAA\n",
	}
	readers := map[string]io.Reader{}
	for k, v := range inputs {
		readers[k] = strings.NewReader(v)
	}
	results, err := svc.AlignBatch(context.Background(), readers, nil)
	if err != nil {
		t.Fatalf("batch failed: %v", err)
	}
	if len(results) != 2 {
		t.Fatalf("want 2 results, got %d", len(results))
	}
	if _, ok := results["job-a"]; !ok {
		t.Fatal("job-a missing")
	}
}

func TestPipelineBatchIsolatesFailures(t *testing.T) {
	svc := New(DefaultConfig())
	readers := map[string]io.Reader{
		"good": strings.NewReader(">a\nACGTACGT\n>b\nACGTACGA\n"),
		"bad":  strings.NewReader(">a\nZZZZ\n"),
	}
	results, err := svc.AlignBatch(context.Background(), readers, nil)
	if err == nil {
		t.Fatal("expected batch error")
	}
	if _, ok := results["good"]; !ok {
		t.Fatal("good job must still be returned")
	}
	if _, ok := results["bad"]; ok {
		t.Fatal("bad job must not appear in results")
	}
}
