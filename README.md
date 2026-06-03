# blast101-cli

Command-line interface and automated test suite extending the BLAST 101
educational protein sequence search application.

**MSc Bioinformatics — University of Edinburgh**  
Bioinformatics Algorithms, 2025–2026

---

## Overview

BLAST 101 is an educational implementation of the BLAST heuristic search
algorithm (Tomlinson, University of Edinburgh, 2025). This repository contains
two components authored as a coursework extension:

- A **command-line interface** that validates inputs, configures parameters,
  and dispatches execution across three search modes from a single entry point
- An **automated test suite** that verifies correctness of the Smith-Waterman
  scoring module and the full BLAST search pipeline

> **Note:** The base BLAST 101 search engine (`blast_101_search.py`,
> `smith_waterman_p.py`, and supporting modules) was written by
> Dr Simon Tomlinson (University of Edinburgh) and is not included here.
> This repository contains only the components authored as part of the
> coursework extension. The full application requires the complete course
> codebase to run.

---

## Repository contents

| File | Description |
|------|-------------|
| `blast_101_cli.py` | CLI entry point: argument parsing, input validation, mode dispatcher |
| `blast_101_tests.py` | Automated unit and integration tests |

---

## Usage

```bash
# BLAST search (default)
python3 blast_101_cli.py -q MNAGWTL -d uniprot_bit2.fasta

# Smith-Waterman exhaustive alignment
python3 blast_101_cli.py -q MNAGWTL -d uniprot_bit2.fasta --mode sw

# Statistical analysis
python3 blast_101_cli.py --mode stats

# Run automated tests
python3 blast_101_cli.py --test
```

### Arguments

| Flag | Required | Description |
|------|----------|-------------|
| `-q` / `--query` | Yes* | Protein query sequence |
| `-d` / `--database` | Yes* | Path to FASTA database file |
| `-m` / `--mode` | No | `blast` (default), `sw`, or `stats` |
| `-t` / `--test` | No | Run test suite and exit |

*Not required for `--mode stats` or `--test`

---

## Input validation

The CLI validates both inputs before any search begins:

**Sequence validation** (`-q`):
- Rejects empty strings
- Detects likely DNA input using protein-exclusive character detection
- Validates each residue against the 20 standard amino acids

**Database validation** (`-d`):
- Confirms the file exists and is non-empty
- Verifies FASTA format (first non-empty line starts with `>`)

---

## Test suite

Six automated tests across two modules, run without manual intervention:

| Test | Type | Ground truth |
|------|------|--------------|
| SW-1: identical sequence (MARIQ) | Unit | BLOSUM62 diagonal sum |
| SW-2: single mismatch (MARIQ vs MARIE) | Unit | BLOSUM62 position scores |
| SW-3: local alignment (QY in MARIQYMARIA) | Unit | BLOSUM62 analytical |
| BLAST-1: hit returned for known query | Integration | Sequence present in DB |
| BLAST-2: top hit is CYP7A1_HUMAN | Integration | Sequence present in DB |
| BLAST-3/4: top hit matches NCBI blastp | Integration | NCBI blastp vs SwissProt (April 2026) |

```bash
python3 blast_101_tests.py
```

Expected output:
```
==================================================
  TOTAL: 6/6 tests passed
==================================================
```

---

## Design notes

**Parameter separation:** Query sequence, database, and mode are passed via
CLI and change per run. Algorithmic parameters (BLOSUM matrix version, word
size, gap penalty) remain in `settings.ini` — mirroring the NCBI BLAST web
interface, where users provide a sequence and database while default parameters
are applied transparently.

**DNA detection:** The standard amino acid and DNA alphabets share four
characters (A, C, G, T), making direct comparison insufficient. Validation
searches for protein-exclusive characters; their absence flags likely DNA input.

**Test ground truth strategy:** SW tests use scores derived analytically from
BLOSUM62 (mathematically exact). BLAST tests compare top-hit identity against
NCBI blastp rather than raw scores, since score values differ between
implementations by design.

---

## Requirements

```
blosum>=1.2
```

```bash
pip install -r requirements.txt
```

---

## Author

Mario Antonio Rodriguez Diaz  
MSc Bioinformatics, University of Edinburgh  
Chevening Scholar 2025–2026  
[linkedin.com/in/mario-rd](https://linkedin.com/in/mario-rd) · [ORCID: 0009-0008-0104-7421](https://orcid.org/0009-0008-0104-7421)
