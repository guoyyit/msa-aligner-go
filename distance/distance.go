package distance

import (
	"errors"
	"math"
	"runtime"
	"sort"
	"sync"

	"github.com/guoyyit/msa-aligner-go/fasta"
)

var errNoSequences = errors.New("distance: no sequences")

// DefaultK is used when Options.K is zero.
const DefaultK = 3

// Alphabet names accepted by Options.
const (
	AlphabetDNA     = fasta.AlphabetDNA
	AlphabetRNA     = fasta.AlphabetRNA
	AlphabetProtein = fasta.AlphabetProtein
)

// Options controls pairwise distance computation.
type Options struct {
	// K is the k-mer length. Zero selects the package default.
	K int
	// Concurrency limits parallel pair computations. Zero selects the default.
	Concurrency int
	// Alphabet controls strand-canonicalization; empty means DNA.
	Alphabet string
}

// Matrix is a symmetric N*N distance matrix stored row-major.
// A failed pair is recorded in Errors and its cell holds NaN.
type Matrix struct {
	N      int
	IDs    []string
	Cells  []float64
	Errors map[[2]int]string
}

// ProgressFunc reports completed/total pair computations.
type ProgressFunc func(done, total int)

// NewMatrix allocates a distance matrix for the given IDs.
// Off-diagonal cells start as NaN; failed pairs remain NaN.
func NewMatrix(ids []string) *Matrix {
	n := len(ids)
	cells := make([]float64, n*n)
	for i := range cells {
		cells[i] = math.NaN()
	}
	for i := 0; i < n; i++ {
		cells[i*n+i] = 0
	}
	return &Matrix{N: n, IDs: ids, Cells: cells, Errors: map[[2]int]string{}}
}

// At returns the distance between leaves i and j.
func (m *Matrix) At(i, j int) float64 { return m.Cells[i*m.N+j] }

// Failed reports whether pair (i,j) failed and why.
func (m *Matrix) Failed(i, j int) (string, bool) {
	if i > j {
		i, j = j, i
	}
	msg, ok := m.Errors[[2]int{i, j}]
	return msg, ok
}

// OrderedErrors returns failed pairs in deterministic index order.
func (m *Matrix) OrderedErrors() []struct {
	I, J int
	Msg  string
} {
	keys := make([][2]int, 0, len(m.Errors))
	for k := range m.Errors {
		keys = append(keys, k)
	}
	sort.Slice(keys, func(a, b int) bool {
		if keys[a][0] != keys[b][0] {
			return keys[a][0] < keys[b][0]
		}
		return keys[a][1] < keys[b][1]
	})
	out := make([]struct {
		I, J int
		Msg  string
	}, len(keys))
	for i, k := range keys {
		out[i] = struct {
			I, J int
			Msg  string
		}{k[0], k[1], m.Errors[k]}
	}
	return out
}

type kmerProfile struct {
	counts map[uint64]int
	total  int
	err    string
}

var dnaComplement = map[byte]byte{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

func buildProfile(seq []byte, k int, alphabet string) kmerProfile {
	counts := make(map[uint64]int)
	if len(seq) < k {
		return kmerProfile{counts: counts, err: "sequence shorter than k"}
	}
	dna := alphabet == "" || alphabet == AlphabetDNA
	for start := 0; start+k <= len(seq); start++ {
		window := seq[start : start+k]
		var key, rc uint64
		usable := true
		for i, c := range window {
			if dna {
				comp, ok := dnaComplement[c]
				if !ok {
					usable = false
					break
				}
				key = key<<2 | dnaBit(c)
				rc = rc | uint64(dnaBit(comp))<<(uint(i)*2)
			} else {
				key = key*256 | uint64(c)
			}
		}
		if !usable {
			continue
		}
		if dna && rc < key {
			key = rc
		}
		counts[key]++
	}
	total := 0
	for _, c := range counts {
		total += c
	}
	p := kmerProfile{counts: counts, total: total}
	if total == 0 {
		p.err = "no usable k-mers"
	}
	return p
}

func dnaBit(c byte) uint64 {
	switch c {
	case 'A':
		return 0
	case 'C':
		return 1
	case 'G':
		return 2
	default:
		return 3
	}
}

func pairDistance(a, b kmerProfile) (float64, string) {
	if a.err != "" {
		return math.NaN(), a.err
	}
	if b.err != "" {
		return math.NaN(), b.err
	}
	common := 0
	for key, ca := range a.counts {
		if cb, ok := b.counts[key]; ok {
			if ca < cb {
				common += ca
			} else {
				common += cb
			}
		}
	}
	denom := a.total
	if b.total < denom {
		denom = b.total
	}
	return 1 - float64(common)/float64(denom), ""
}

// Compute calculates every pairwise distance with bounded concurrency.
// A single failing pair never aborts the remaining pairs.
func Compute(seqs []fasta.Sequence, opts Options, progress ProgressFunc) (*Matrix, error) {
	if len(seqs) == 0 {
		return nil, errNoSequences
	}
	k := opts.K
	if k <= 0 {
		k = DefaultK
	}
	workers := opts.Concurrency
	if workers <= 0 {
		workers = runtime.NumCPU()
	}

	ids := make([]string, len(seqs))
	profiles := make([]kmerProfile, len(seqs))
	for i, s := range seqs {
		ids[i] = s.ID
		profiles[i] = buildProfile(s.Data, k, opts.Alphabet)
	}
	matrix := NewMatrix(ids)

	type pair struct{ i, j int }
	total := len(seqs) * (len(seqs) - 1) / 2
	type result struct {
		i, j int
		d    float64
		msg  string
	}
	tasks := make(chan pair)
	results := make(chan result, total)
	var wg sync.WaitGroup
	for w := 0; w < workers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for t := range tasks {
				d, msg := pairDistanceSafe(profiles[t.i], profiles[t.j])
				results <- result{t.i, t.j, d, msg}
			}
		}()
	}
	go func() {
		for i := 0; i < len(seqs); i++ {
			for j := i + 1; j < len(seqs); j++ {
				tasks <- pair{i, j}
			}
		}
		close(tasks)
	}()
	doneCh := make(chan struct{})
	go func() {
		wg.Wait()
		close(doneCh)
	}()

	done := 0
collect:
	for {
		select {
		case r := <-results:
			matrix.Cells[r.i*matrix.N+r.j] = r.d
			matrix.Cells[r.j*matrix.N+r.i] = r.d
			if r.msg != "" {
				matrix.Errors[[2]int{r.i, r.j}] = r.msg
			}
			done++
			if progress != nil {
				progress(done, total)
			}
		case <-doneCh:
			// Drain any results delivered concurrently with completion.
			for {
				select {
				case r := <-results:
					matrix.Cells[r.i*matrix.N+r.j] = r.d
					matrix.Cells[r.j*matrix.N+r.i] = r.d
					if r.msg != "" {
						matrix.Errors[[2]int{r.i, r.j}] = r.msg
					}
					done++
					if progress != nil {
						progress(done, total)
					}
				default:
					break collect
				}
			}
		}
	}
	return matrix, nil
}

func pairDistanceSafe(a, b kmerProfile) (d float64, msg string) {
	defer func() {
		if r := recover(); r != nil {
			d = math.NaN()
			msg = "pair computation panicked"
		}
	}()
	return pairDistance(a, b)
}

// PairDistance computes the distance of one sequence pair with the same
// semantics Compute applies (k-mer profile distance, NaN + reason on
// unusable input). It exists for resuming individual missing pairs.
func PairDistance(a, b []byte, opts Options) (float64, string) {
	k := opts.K
	if k <= 0 {
		k = DefaultK
	}
	pa := buildProfile(a, k, opts.Alphabet)
	pb := buildProfile(b, k, opts.Alphabet)
	return pairDistanceSafe(pa, pb)
}
