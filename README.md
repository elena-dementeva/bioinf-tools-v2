# Bioinf Tools

**Bioinf Tools** is a Python toolkit designed for DNA/RNA sequence manipulations and FASTQ sequence filtering based on GC content, sequence length and quality threshold. The toolkit provides two main functionalities:
1. **run_dna_rna_tools**: Manipulate DNA/RNA sequences (transcribe, reverse, complement and generate reverse complements).
2. **filter_fastq**: Filter FASTQ sequences based on defined thresholds for GC content, length, and quality.
3. **convert_multiline_fasta_to_oneline**: Convert multi-line FASTA sequences into single-line sequences.
4. **parse_blast_output**: Extract best match descriptions from BLAST output files.
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
- `input_fastq`: input FASTQ file.
- `output_fastq`: filtered sequences (saved in the filtered folder).
- `gc_bounds`: The GC content bounds for filtering. Default is (0, 100).
- `length_bounds`: The sequence length bounds for filtering. Default is (0, 2**32).
- `quality_threshold`: The threshold for filtering by sequence quality. Default is 0.

#### Returns:
- Filtered sequences that meet the filtering criteria are written to the specified output_fastq file in the filtered folder.

---

### 3. `convert_multiline_fasta_to_oneline`

Converts multi-line FASTA sequences to single-line sequences, saving each sequence in a single line in a new FASTA file.

#### Arguments:
- `input_fasta`: input FASTA file.
- `output_fasta`: output FASTA file.

#### Returns:
- Writes the reformatted sequences to a specified FASTA file.

---

### 4. `parse_blast_output`

Extracts the best match descriptions from a BLAST output file and saves them in alphabetical order in a new file.

#### Arguments:
- `input_file`: input BLAST file.
- `output_file`: extracted descriptions.

#### Returns:
- Writes the extracted and sorted descriptions to a output file.

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

# Filter sequences in the example_fastq.fastq file
filter_fastq("example_fastq.fastq", "filtered_output.fastq", gc_bounds=(40, 60), length_bounds=(6, 10), quality_threshold=35)
```
This filters sequences from the specified file based on the given thresholds and saves them to filtered folder.

### Example 3: Convert Multi-line FASTA to Single-line

```python
from bio_files_processor import convert_multiline_fasta_to_oneline

# Convert multi-line FASTA to single-line FASTA
convert_multiline_fasta_to_oneline("input.fasta", "output_oneline.fasta")
```
### Example 4: Parse BLAST Output

```python
from bio_files_processor import parse_blast_output

# Parse BLAST output and save best match descriptions
parse_blast_output("example_blast_results.txt", "parsed_results.txt")
```

## Contact
For any questions or feedback please contact:

Elena Dementeva
Email: dementevayelena@gmail.com

