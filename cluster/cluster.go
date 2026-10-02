package cluster

import (
	"errors"
	"fmt"
	"math"

	"github.com/guoyyit/msa-aligner-go/distance"
)

var (
	errNoLeaves   = errors.New("cluster: distance matrix has no leaves")
	errAllFailed  = errors.New("cluster: all pairwise distances failed")
	errMatrixSize = errors.New("cluster: inconsistent matrix size")
)

// MergeStep describes one deterministic agglomerative merge.
type MergeStep struct {
	Step     int
	Left     int
	Right    int
	Distance float64
	NewIndex int
}

// Tree is the fully ordered merge plan.
type Tree struct {
	N     int
	IDs   []string
	Steps []MergeStep
	// Members[node] lists the original leaf indices contained in node.
	// Nodes 0..N-1 are leaves; merged nodes are appended in step order.
	Members [][]int
}

// Build derives the unique merge order from the distance matrix.
//
// Linkage is average linkage over all *finite* leaf-pair distances:
// clusters are weighted by the number of finite pairs they contribute.
// Pairs with no finite distance (NaN) are considered unavailable and are
// only selected after every finite candidate has been exhausted, in
// ascending node-index order. Ties on distance are broken by the pair
// (smaller node, larger node), so identical input always yields the
// identical merge order on every machine.
func Build(m *distance.Matrix) (*Tree, error) {
	if m == nil || m.N == 0 {
		return nil, errNoLeaves
	}
	if len(m.Cells) != m.N*m.N {
		return nil, errMatrixSize
	}
	n := m.N
	active := make(map[int]bool, n)
	type agg struct {
		sum float64
		cnt int
	}
	acc := make(map[[2]int]*agg)
	for i := 0; i < n; i++ {
		active[i] = true
	}
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			d := m.At(i, j)
			if !math.IsNaN(d) {
				acc[[2]int{i, j}] = &agg{sum: d, cnt: 1}
			}
		}
	}
	if n > 1 && len(acc) == 0 {
		return nil, errAllFailed
	}

	members := make([][]int, 2*n-1)
	for i := 0; i < n; i++ {
		members[i] = []int{i}
	}
	tree := &Tree{N: n, IDs: append([]string(nil), m.IDs...), Members: members}

	nextNode := n
	for len(active) > 1 {
		bestA, bestB, bestD, bestFinite := -1, -1, math.NaN(), false
		for a := range active {
			for b := range active {
				if a >= b {
					continue
				}
				entry, finite := acc[pairKey(a, b)]
				if !finite || entry.cnt == 0 {
					continue
				}
				d := entry.sum / float64(entry.cnt)
				if !bestFinite || d < bestD || (d == bestD && pairBefore(a, b, bestA, bestB)) {
					bestA, bestB, bestD, bestFinite = a, b, d, true
				}
			}
		}
		finiteChosen := true
		if !bestFinite {
			// No finite candidate remains: deterministically pick the
			// smallest-indexed active pair; its recorded distance is NaN.
			finiteChosen = false
			for a := range active {
				for b := range active {
					if a < b && (bestA == -1 || pairBefore(a, b, bestA, bestB)) {
						bestA, bestB = a, b
					}
				}
			}
		}

		left, right := bestA, bestB
		if left > right {
			left, right = right, left
		}
		newNode := nextNode
		nextNode++
		step := MergeStep{
			Step:     len(tree.Steps),
			Left:     left,
			Right:    right,
			Distance: bestD,
			NewIndex: newNode,
		}
		if !finiteChosen {
			step.Distance = math.NaN()
		}
		tree.Steps = append(tree.Steps, step)
		tree.Members[newNode] = append(append([]int(nil), members[left]...), members[right]...)

		delete(active, left)
		delete(active, right)
		active[newNode] = true
		for x := range active {
			if x == newNode {
				continue
			}
			s, c := 0.0, 0
			if e := acc[pairKey(left, x)]; e != nil {
				s += e.sum
				c += e.cnt
			}
			if e := acc[pairKey(right, x)]; e != nil {
				s += e.sum
				c += e.cnt
			}
			if c == 0 {
				continue
			}
			acc[pairKey(newNode, x)] = &agg{sum: s, cnt: c}
		}
	}

	return tree, nil
}

func pairKey(a, b int) [2]int {
	if a > b {
		a, b = b, a
	}
	return [2]int{a, b}
}

func pairBefore(a, b, c, d int) bool {
	if a != c {
		return a < c
	}
	return b < d
}

func (t *Tree) String() string {
	return fmt.Sprintf("cluster.Tree(N=%d, merges=%d)", t.N, len(t.Steps))
}
