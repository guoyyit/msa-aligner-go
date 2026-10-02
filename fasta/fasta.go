package fasta

import (
	"bufio"
	"fmt"
	"io"
	"strings"
)

// Sequence is one parsed FASTA record.
type Sequence struct {
	ID          string
	Description string
	Data        []byte
}

// ErrorKind categorizes FASTA parse failures so callers can distinguish them.
type ErrorKind string

const (
	ErrKindMissingID        ErrorKind = "missing_id"
	ErrKindEmptySequence    ErrorKind = "empty_sequence"
	ErrKindIllegalCharacter ErrorKind = "illegal_character"
	ErrKindDuplicateID      ErrorKind = "duplicate_id"
	ErrKindMalformed        ErrorKind = "malformed"
)

// Supported alphabet names used by Parser.
const (
	AlphabetDNA     = "dna"
	AlphabetRNA     = "rna"
	AlphabetProtein = "protein"
)

var alphabets = map[string]map[byte]byte{
	AlphabetDNA:     buildAlphabet("ACGTN-", "T"),
	AlphabetRNA:     buildAlphabet("ACGUN-", "U"),
	AlphabetProtein: buildAlphabet("ACDEFGHIKLMNPQRSTVWYBXZ*-", ""),
}

func buildAlphabet(chars, rnaU string) map[byte]byte {
	m := make(map[byte]byte, len(chars))
	for i := 0; i < len(chars); i++ {
		c := chars[i]
		m[c] = c
		if c >= 'A' && c <= 'Z' {
			m[c+('a'-'A')] = c
		}
	}
	if rnaU == "T" { // DNA input normalizes U to T.
		m['U'] = 'T'
		m['u'] = 'T'
	}
	if rnaU == "U" { // RNA input normalizes T to U.
		m['T'] = 'U'
		m['t'] = 'U'
	}
	return m
}

// ParseError describes a single rejected FASTA input with a stable Kind.
type ParseError struct {
	Kind   ErrorKind
	ID     string
	Line   int
	Column int
	Detail string
}

func (e *ParseError) Error() string {
	where := ""
	if e.ID != "" {
		where = " record " + e.ID
	}
	if e.Line > 0 {
		where += fmt.Sprintf(" at line %d", e.Line)
		if e.Column > 0 {
			where += fmt.Sprintf(" column %d", e.Column)
		}
	}
	detail := ""
	if e.Detail != "" {
		detail = ": " + e.Detail
	}
	return fmt.Sprintf("fasta %s%s%s", e.Kind, where, detail)
}

// Errors collects every independently detectable parse failure.
type Errors []*ParseError

func (es Errors) Error() string {
	parts := make([]string, len(es))
	for i, e := range es {
		parts[i] = e.Error()
	}
	return strings.Join(parts, "; ")
}

// Parser parses FASTA streams against one alphabet.
type Parser struct {
	// Alphabet is one of the Alphabet* constants; empty means AlphabetDNA.
	Alphabet string
}

// Parse reads FASTA records using the DNA alphabet.
func Parse(r io.Reader) ([]Sequence, error) {
	return Parser{}.Parse(r)
}

// Parse reads and fully validates the FASTA stream.
// Failures are returned as Errors containing typed ParseError values.
func (p Parser) Parse(r io.Reader) ([]Sequence, error) {
	alphabet := p.Alphabet
	if alphabet == "" {
		alphabet = AlphabetDNA
	}
	legal, ok := alphabets[alphabet]
	if !ok {
		return nil, &ParseError{Kind: ErrKindMalformed, Detail: "unknown alphabet: " + alphabet}
	}

	var errs Errors
	var seqs []Sequence
	seen := make(map[string]int)

	var curID, curDesc string
	var curData []byte
	headerLine := 0
	haveRecord := false

	finish := func() {
		if !haveRecord {
			return
		}
		if len(curData) == 0 {
			errs = append(errs, &ParseError{Kind: ErrKindEmptySequence, ID: curID, Line: headerLine})
			return
		}
		if prev, dup := seen[curID]; dup {
			errs = append(errs, &ParseError{
				Kind:   ErrKindDuplicateID,
				ID:     curID,
				Line:   headerLine,
				Detail: fmt.Sprintf("first seen at record #%d", prev+1),
			})
			return
		}
		seen[curID] = len(seqs)
		seqs = append(seqs, Sequence{ID: curID, Description: curDesc, Data: curData})
	}

	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
	lineNo := 0
	for scanner.Scan() {
		lineNo++
		line := scanner.Text()
		trimmed := strings.TrimRight(line, "\r\n")
		if strings.TrimSpace(trimmed) == "" {
			continue
		}
		if strings.HasPrefix(trimmed, ">") {
			finish()
			haveRecord = true
			headerLine = lineNo
			header := strings.TrimSpace(strings.TrimPrefix(trimmed, ">"))
			if header == "" {
				curID, curDesc, curData = "", "", nil
				errs = append(errs, &ParseError{Kind: ErrKindMissingID, Line: lineNo})
				continue
			}
			// First whitespace-delimited token is the ID; remainder is description.
			fields := strings.Fields(header)
			curID = fields[0]
			if len(fields) > 1 {
				curDesc = strings.TrimSpace(strings.TrimPrefix(header, fields[0]))
			} else {
				curDesc = ""
			}
			curData = nil
			continue
		}
		if !haveRecord {
			errs = append(errs, &ParseError{Kind: ErrKindMalformed, Line: lineNo, Detail: "sequence data before first '>' header"})
			continue
		}
		for col, raw := range []byte(trimmed) {
			if raw == ' ' || raw == '\t' {
				continue
			}
			canonical, legalChar := legal[raw]
			if !legalChar {
				errs = append(errs, &ParseError{
					Kind:   ErrKindIllegalCharacter,
					ID:     curID,
					Line:   lineNo,
					Column: col + 1,
					Detail: fmt.Sprintf("character %q not in %s alphabet", raw, alphabet),
				})
				continue
			}
			curData = append(curData, canonical)
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, &ParseError{Kind: ErrKindMalformed, Detail: err.Error()}
	}
	finish()

	if len(seqs) == 0 && len(errs) == 0 {
		errs = append(errs, &ParseError{Kind: ErrKindMalformed, Detail: "input contains no FASTA records"})
	}
	if len(errs) > 0 {
		return seqs, errs
	}
	return seqs, nil
}
