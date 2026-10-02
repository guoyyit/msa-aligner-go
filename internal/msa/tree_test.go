package msa

import (
	"errors"
	"testing"
)

func TestBuildMergeOrderUPGMA(t *testing.T) {
	// a,b close; c,d close; the two groups equidistant.
	dm := &DistanceMatrix{
		N: 4,
		Dist: [][]float64{
			{0, 0.1, 0.5, 0.5},
			{0.1, 0, 0.5, 0.5},
			{0.5, 0.5, 0, 0.2},
			{0.5, 0.5, 0.2, 0},
		},
	}
	order, err := BuildMergeOrder(dm)
	if err != nil {
		t.Fatal(err)
	}
	want := []MergeEvent{{0, 1}, {2, 3}, {4, 5}}
	if len(order.Events) != len(want) {
		t.Fatalf("events = %v, want %v", order.Events, want)
	}
	for i := range want {
		if order.Events[i] != want[i] {
			t.Errorf("event %d = %+v, want %+v", i, order.Events[i], want[i])
		}
	}
}

func TestBuildMergeOrderDeterministic(t *testing.T) {
	dm := &DistanceMatrix{
		N: 4,
		Dist: [][]float64{
			{0, 0.3, 0.2, 0.4},
			{0.3, 0, 0.5, 0.1},
			{0.2, 0.5, 0, 0.6},
			{0.4, 0.1, 0.6, 0},
		},
	}
	o1, _ := BuildMergeOrder(dm)
	o2, _ := BuildMergeOrder(dm)
	if len(o1.Events) != len(o2.Events) {
		t.Fatalf("lengths differ")
	}
	for i := range o1.Events {
		if o1.Events[i] != o2.Events[i] {
			t.Errorf("event %d differs: %+v vs %+v", i, o1.Events[i], o2.Events[i])
		}
	}
}

func TestBuildMergeOrderTieBreak(t *testing.T) {
	// All pairwise distances equal -> ties broken lexicographically.
	dm := &DistanceMatrix{
		N: 3,
		Dist: [][]float64{
			{0, 0.5, 0.5},
			{0.5, 0, 0.5},
			{0.5, 0.5, 0},
		},
	}
	order, err := BuildMergeOrder(dm)
	if err != nil {
		t.Fatal(err)
	}
	// First merge must be (0,1) lexicographically.
	if order.Events[0] != (MergeEvent{0, 1}) {
		t.Errorf("first merge = %+v, want {0 1}", order.Events[0])
	}
}

func TestBuildMergeOrderSingle(t *testing.T) {
	dm := &DistanceMatrix{N: 1, Dist: [][]float64{{0}}}
	order, err := BuildMergeOrder(dm)
	if err != nil {
		t.Fatal(err)
	}
	if len(order.Events) != 0 {
		t.Errorf("events = %v, want empty", order.Events)
	}
}

func TestBuildMergeOrderFailedPair(t *testing.T) {
	dm := &DistanceMatrix{
		N: 3,
		Dist: [][]float64{
			{0, -1, 0.5},
			{-1, 0, 0.5},
			{0.5, 0.5, 0},
		},
	}
	_, err := BuildMergeOrder(dm)
	if !errors.Is(err, ErrIncompleteDistance) {
		t.Errorf("err = %v, want ErrIncompleteDistance", err)
	}
}
