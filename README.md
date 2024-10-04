# Bioinf Tools

**Bioinf Tools** is a Python toolkit designed for DNA/RNA sequence manipulations and FASTQ sequence filtering based on GC content, sequence length and quality threshold. The toolkit provides two main functionalities:
1. **run_dna_rna_tools**: Manipulate DNA/RNA sequences (transcribe, reverse, complement and generate reverse complements).
2. **filter_fastq**: Filter fastq sequences based on defined thresholds for GC content, length, and quality.

---

## Authors:

- **Software**: Elena Dementeva (dementevayelena@gmail.com)

---

## Content

- [Installation](#installation)
- [Functions Overview](#functions-overview)
- [Usage Examples](#usage-examples)
- [Contact](#contact)

---

## Installation

To use **Bioinf Tools**, simply clone the repository to your local machine:

```bash
git clone https://github.com/elena-dementeva/bioinf-tools-v2
cd bioinf-tools
```
Ensure you have Python 3.x installed. No additional external libraries are required.

## Functions Overview

### 1. `run_dna_rna_tools`

This function allows users to perform DNA/RNA sequence manipulations, including:
- **transcribe**: Convert DNA to RNA.
- **reverse**: Reverse a sequence.
- **complement**: Generate the complement of a sequence.
- **reverse_complement**: Generate the reverse complement.

#### Arguments:
- `*args`: The sequences followed by the procedure type (`"transcribe"`, `"reverse"`, `"complement"`, `"reverse_complement"`).

#### Returns:
- A single manipulated sequence (if one is provided) or a list of sequences.

---

### 2. `filter_fastq`

This function filters fastq sequences based on the following criteria: GC content, sequence length and quality.

#### Arguments:
- `seqs`: is a dictionary where:
  - The key is the sequence name.
  - The value is a tuple with the sequence and the quality scores.
- `gc_bounds`: The GC content bounds for filtering. Default is (0, 100).
- `length_bounds`: The sequence length bounds for filtering. Default is (0, 2**32).
- `quality_threshold`: The threshold for filtering by sequence quality. Default is 0.

#### Returns:
- A dictionary with sequences that meet the filtering criteria.

---

## Usage Examples

### Example 1: Manipulate DNA/RNA Sequences

```python
from bioinf_tools import run_dna_rna_tools

# Transcribe DNA to RNA
result = run_dna_rna_tools("ATCG", "GCTA", "transcribe")
print(result)
```
This transcribes both provided sequences.

### Example 2: Filter Fastq Sequences

```python
from bioinf_tools import filter_fastq

# Example fastq sequences
fastq_data = {
    'seq1': ('ATCGATCG', 'IIIIIIII'),
    'seq2': ('GGCCGGCC', 'HHHHHHHH'),
    'seq3': ('TTAAttaa', 'IIIIHHHH'),
}

# Filter sequences by GC content and quality
filtered = filter_fastq(fastq_data, gc_bounds=(40, 60), length_bounds=(6, 10), quality_threshold=35)
print(filtered)
```
This filters sequences from the provided fastq_data based on the given thresholds.

## Contact
For any questions or feedback please contact:

Elena Dementeva
Email: dementevayelena@gmail.com

