package msa

import (
	"fmt"
	"sort"
	"sync"
)

// DistanceMatrix is an n x n symmetric matrix of pairwise distances in [0,1].
type DistanceMatrix struct {
	N    int
	Dist [][]float64
}

// ComputeDistances computes all pairwise distances concurrently.
// A failure on one pair is isolated and reported in the returned error list;
// it never panics and never blocks other pairs. Failed pairs are set to -1.
func ComputeDistances(seqs []Sequence, k int, concurrency int, prog ProgressFunc) (*DistanceMatrix, []error) {
	n := len(seqs)
	dm := &DistanceMatrix{N: n, Dist: make([][]float64, n)}
	for i := range dm.Dist {
		dm.Dist[i] = make([]float64, n)
	}
	if concurrency < 1 {
		concurrency = 1
	}

	type job struct{ i, j int }
	jobs := make(chan job)
	var wg sync.WaitGroup
	var mu sync.Mutex
	var errs []error
	var done int64
	total := n * (n - 1) / 2

	worker := func() {
		defer wg.Done()
		for jb := range jobs {
			func() {
				defer func() {
					if r := recover(); r != nil {
						mu.Lock()
						errs = append(errs, fmt.Errorf("distance pair (%d,%d) panic: %v", jb.i, jb.j, r))
						dm.Dist[jb.i][jb.j] = -1
						dm.Dist[jb.j][jb.i] = -1
						done++
						mu.Unlock()
					}
				}()
				d, err := PairDistance(seqs[jb.i].Seq, seqs[jb.j].Seq, k)
				mu.Lock()
				if err != nil {
					errs = append(errs, fmt.Errorf("distance pair (%d,%d): %w", jb.i, jb.j, err))
					dm.Dist[jb.i][jb.j] = -1
					dm.Dist[jb.j][jb.i] = -1
				} else {
					dm.Dist[jb.i][jb.j] = d
					dm.Dist[jb.j][jb.i] = d
				}
				done++
				cur := int(done)
				mu.Unlock()
				if prog != nil {
					prog(Progress{
						Phase:   "distance",
						Current: cur,
						Total:   total,
						Message: fmt.Sprintf("computed %d/%d pairs", done, total),
					})
				}
			}()
		}
	}

	wg.Add(concurrency)
	for w := 0; w < concurrency; w++ {
		go worker()
	}
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			jobs <- job{i, j}
		}
	}
	close(jobs)
	wg.Wait()
	return dm, errs
}

// kmerCounts returns overlapping k-mer frequencies.
func kmerCounts(s string, k int) map[string]int {
	counts := make(map[string]int)
	if len(s) < k || k < 1 {
		return counts
	}
	for i := 0; i <= len(s)-k; i++ {
		counts[s[i:i+k]]++
	}
	return counts
}

// PairDistance computes the k-mer (Dice) distance between two sequences.
// For sequences shorter than k it degrades to residue-composition (k=1).
func PairDistance(a, b string, k int) (float64, error) {
	if k < 1 {
		k = 1
	}
	ka := kmerCounts(a, k)
	kb := kmerCounts(b, k)
	if len(ka) == 0 || len(kb) == 0 {
		ka = kmerCounts(a, 1)
		kb = kmerCounts(b, 1)
	}

	var totalA, totalB float64
	keys := make([]string, 0, len(ka))
	for key, c := range ka {
		keys = append(keys, key)
		totalA += float64(c)
	}
	for _, c := range kb {
		totalB += float64(c)
	}
	sort.Strings(keys)

	var shared float64
	for _, key := range keys {
		shared += float64(min(ka[key], kb[key]))
	}
	if totalA+totalB == 0 {
		return 0, fmt.Errorf("%w: cannot distance empty sequences", ErrIncompleteDistance)
	}
	dice := 2 * shared / (totalA + totalB)
	if dice > 1 {
		dice = 1
	}
	return 1 - dice, nil
}
