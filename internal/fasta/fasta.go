package fasta

import (
	"errors"
	"fmt"
	"strings"
	"unicode"
)

// Sequence is a single FASTA record.
type Sequence struct {
	ID   string // identifier (first token after '>')
	Desc string // free-form description after the identifier
	Seq  string // residues (uppercase); gaps only in aligned input
}

// Validation errors with distinguishable failure reasons.
var (
	ErrEmptySequence    = errors.New("fasta: empty sequence")
	ErrInvalidCharacter  = errors.New("fasta: invalid character in sequence")
	ErrDuplicateID       = errors.New("fasta: duplicate sequence identifier")
	ErrEmptyInput        = errors.New("fasta: no sequences found")
	ErrEmptyHeader       = errors.New("fasta: empty header")
	ErrSeqWithoutHeader  = errors.New("fasta: sequence data before any header")
	ErrRaggedAlignment   = errors.New("fasta: aligned rows have unequal lengths")
)

// validResidues is the union of IUPAC DNA/RNA and protein residue codes.
var validResidues = map[rune]bool{
	'A': true, 'B': true, 'C': true, 'D': true, 'E': true, 'F': true,
	'G': true, 'H': true, 'I': true, 'K': true, 'L': true, 'M': true,
	'N': true, 'P': true, 'Q': true, 'R': true, 'S': true, 'T': true,
	'U': true, 'V': true, 'W': true, 'X': true, 'Y': true, 'Z': true,
	'*': true,
}

// Parse reads FASTA-formatted raw sequences and returns validated sequences.
// Empty sequences, invalid characters, and duplicate identifiers are rejected
// with distinguishable errors. Gap characters are rejected in raw input.
func Parse(data []byte) ([]Sequence, error) {
	return parse(data, false)
}

// ParseAlignment reads an aligned FASTA (gap characters '-' and '.' allowed).
// All rows must have equal length.
func ParseAlignment(data []byte) ([]Sequence, error) {
	seqs, err := parse(data, true)
	if err != nil {
		return nil, err
	}
	if len(seqs) > 0 {
		want := len(seqs[0].Seq)
		for i, s := range seqs {
			if len(s.Seq) != want {
				return nil, fmt.Errorf("%w: row %d (%s) length %d != %d", ErrRaggedAlignment, i, s.ID, len(s.Seq), want)
			}
		}
	}
	return seqs, nil
}

func parse(data []byte, allowGaps bool) ([]Sequence, error) {
	var seqs []Sequence
	var cur *Sequence
	seen := make(map[string]bool)

	flush := func() error {
		if cur == nil {
			return nil
		}
		if len(cur.Seq) == 0 {
			return fmt.Errorf("%w: %q", ErrEmptySequence, cur.ID)
		}
		if err := validateChars(cur, allowGaps); err != nil {
			return err
		}
		if seen[cur.ID] {
			return fmt.Errorf("%w: %q", ErrDuplicateID, cur.ID)
		}
		seen[cur.ID] = true
		seqs = append(seqs, *cur)
		cur = nil
		return nil
	}

	for _, raw := range strings.Split(string(data), "\n") {
		line := strings.TrimRight(raw, "\r")
		if strings.HasPrefix(line, ">") {
			if err := flush(); err != nil {
				return nil, err
			}
			header := line[1:]
			fields := strings.Fields(header)
			if len(fields) == 0 {
				return nil, fmt.Errorf("%w: %q", ErrEmptyHeader, line)
			}
			id := fields[0]
			idx := strings.Index(header, id)
			desc := ""
			if idx >= 0 {
				desc = strings.TrimSpace(header[idx+len(id):])
			}
			cur = &Sequence{ID: id, Desc: desc}
			continue
		}

		stripped := strings.TrimSpace(line)
		if stripped == "" {
			continue
		}
		if cur == nil {
			return nil, fmt.Errorf("%w: %q", ErrSeqWithoutHeader, stripped)
		}
		cur.Seq += strings.ToUpper(strings.Map(func(r rune) rune {
			if unicode.IsSpace(r) {
				return -1
			}
			return r
		}, line))
	}

	if err := flush(); err != nil {
		return nil, err
	}
	if len(seqs) == 0 {
		return nil, ErrEmptyInput
	}
	return seqs, nil
}

func validateChars(s *Sequence, allowGaps bool) error {
	for i, r := range s.Seq {
		if validResidues[r] {
			continue
		}
		if allowGaps && (r == '-' || r == '.') {
			continue
		}
		return fmt.Errorf("%w: id=%q pos=%d char=%q", ErrInvalidCharacter, s.ID, i+1, r)
	}
	return nil
}
