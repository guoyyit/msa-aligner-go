package msa

import (
	"fmt"
	"math"
)

// MergeEvent records one UPGMA cluster merge.
type MergeEvent struct {
	Left  int // cluster id
	Right int // cluster id
}

// MergeOrder is the deterministic sequence of cluster merges.
// The k-th event (0-indexed) creates cluster id N+k.
type MergeOrder struct {
	N      int
	Events []MergeEvent // length N-1
}

// BuildMergeOrder constructs a deterministic merge order from a distance matrix
// using UPGMA (unweighted arithmetic mean). Ties are broken by cluster index,
// so identical inputs always produce an identical merge order.
func BuildMergeOrder(dm *DistanceMatrix) (*MergeOrder, error) {
	n := dm.N
	if n < 1 {
		return nil, fmt.Errorf("%w: empty distance matrix", ErrIncompleteDistance)
	}
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if dm.Dist[i][j] < 0 {
				return nil, fmt.Errorf("%w: failed pair (%d,%d)", ErrIncompleteDistance, i, j)
			}
		}
	}
	if n == 1 {
		return &MergeOrder{N: 1, Events: nil}, nil
	}

	size := 2*n - 1
	dist := make([][]float64, size)
	for i := range dist {
		dist[i] = make([]float64, size)
	}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			dist[i][j] = dm.Dist[i][j]
		}
	}
	active := make([]bool, size)
	for i := 0; i < n; i++ {
		active[i] = true
	}

	const eps = 1e-12
	events := make([]MergeEvent, 0, n-1)
	nextID := n
	for step := 0; step < n-1; step++ {
		bestI, bestJ := -1, -1
		bestD := math.Inf(1)
		for i := 0; i < size; i++ {
			if !active[i] {
				continue
			}
			for j := i + 1; j < size; j++ {
				if !active[j] {
					continue
				}
				d := dist[i][j]
				if d < bestD-eps {
					bestD, bestI, bestJ = d, i, j
				} else if math.Abs(d-bestD) <= eps {
					// Tie: break by lexicographically smallest (i,j).
					if bestI < 0 || i < bestI || (i == bestI && j < bestJ) {
						bestI, bestJ = i, j
					}
				}
			}
		}
		if bestI < 0 {
			return nil, fmt.Errorf("%w: no pair to merge at step %d", ErrIncompleteDistance, step)
		}
		events = append(events, MergeEvent{Left: bestI, Right: bestJ})
		newID := nextID
		nextID++
		for k := 0; k < size; k++ {
			if !active[k] || k == bestI || k == bestJ {
				continue
			}
			d := (dist[bestI][k] + dist[bestJ][k]) / 2
			dist[newID][k] = d
			dist[k][newID] = d
		}
		active[bestI] = false
		active[bestJ] = false
		active[newID] = true
	}
	return &MergeOrder{N: n, Events: events}, nil
}
