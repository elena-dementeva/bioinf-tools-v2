# Bioinf Tools

**Bioinf Tools** is a Python toolkit designed for DNA/RNA sequence manipulations and FASTQ sequence filtering based on GC content, sequence length, and quality threshold.  
The toolkit provides two main functionalities:

1. **BiologicalSequence classes**: Abstract classes for handling DNA, RNA, and protein sequences.
2. **filter_fastq**: Filters FASTQ sequences based on defined thresholds for GC content, length, and quality.

---

## Authors:

- **Software**: Elena Dementeva (dementevayelena@gmail.com)

---

## Content

- [Installation](#installation)
- [Class Overview](#class-overview)
- [Usage Examples](#usage-examples)
- [Contact](#contact)

---

## Installation

To use **Bioinf Tools**, simply clone the repository to your local machine:

```bash
git clone https://github.com/elena-dementeva/bioinf-tools-v2
cd bioinf-tools-v2
```
Ensure you have Python 3.x installed and install the required dependencies:
```bash
pip install -r requirements.txt
```

---
## Class Overview

### 1. **`BiologicalSequence`**
An **abstract class** that defines a common interface for working with biological sequences. It provides:

- **`__len__()`** – Returns the length of the sequence.
- **`__getitem__()`** – Enables indexing and slicing.
- **`__str__()` / `__repr__()`** – Returns a string representation.
- **`check_alphabet()`** – Ensures the sequence contains valid characters.
- **`to_oneline_fasta(output_fasta)`** – Saves the sequence in a single-line FASTA format.

---

### 2. **`NucleicAcidSequence`**
A **base class** for handling nucleic acid sequences (DNA and RNA). Implements:

- **`complement()`** – Returns the complementary sequence.
- **`reverse()`** – Returns the reversed sequence.
- **`reverse_complement()`** – Returns the reverse complement of the sequence.

---

### 3. **`DNASequence`**
A **class for handling DNA sequences**. Inherits from **`NucleicAcidSequence`** and adds:

- **`transcribe()`** – Converts DNA to RNA.

---

### 4. **`RNASequence`**
A **class for handling RNA sequences**. Inherits from **`NucleicAcidSequence`**.

---

### 5. **`AminoAcidSequence`**
A **class for handling protein sequences**. Implements:

- **`count_hydrophobic_residues()`** – Counts the number of hydrophobic amino acids.

---

### 6. **`BlastParser`**
A **class for parsing BLAST output files**. Implements:

- **`parse_blast(input_file, output_file)`** – Extracts and sorts the best match descriptions from a BLAST output.

---

### 7. **`filter_fastq`**
A **function** that filters FASTQ sequences based on GC content, sequence length, and quality.

#### Arguments:
- **`input_fastq`** – **(str)** Input FASTQ file.
- **`output_fastq`** – **(str)** Output file for filtered sequences.
- **`gc_bounds`** – **(tuple, optional)** The GC content bounds *(default: `(0, 100)`).*
- **`length_bounds`** – **(tuple, optional)** The sequence length bounds *(default: `(0, 2**32)`).*
- **`quality_threshold`** – **(float, optional)** Minimum quality threshold *(default: `0`).*

---

## Usage Examples

### Example 1: Manipulate DNA/RNA Sequences

```python
from bioinf_tools import DNASequence, RNASequence

# Create a DNA sequence
dna = DNASequence("ATGC")

# Get complement
print(dna.complement())  # Output: TACG

# Get reverse complement
print(dna.reverse_complement())  # Output: GCAT

# Transcribe DNA to RNA
rna = dna.transcribe()
print(rna)  # Output: AUGC
```

### Example 2: Filter Fastq Sequences

```python
from bioinf_tools import filter_fastq

# Filter sequences in example.fastq
filter_fastq(
    "example.fastq", 
    "filtered.fastq", 
    gc_bounds=(40, 60), 
    length_bounds=(50, 150), 
    quality_threshold=30
)
```
This filters sequences from the specified file based on the given thresholds and saves them to an output file.

### Example 3: Parse BLAST Output

```python
from bioinf_tools import BlastParser

# Extract best matches from a BLAST output
BlastParser.parse_blast("blast_output.txt", "parsed_results.txt")
```
---
## Command-Line Interface (CLI)

The FastQ Filter Tool now features a command-line interface that allows you to filter FASTQ files directly from the terminal. This CLI leverages Python’s `argparse` module to accept various parameters:

- **--input**: Path to the input FASTQ file.
- **--output**: Path to the output FASTQ file.
- **--gc**: GC content bounds. Provide either one value (upper bound) or two values (min and max).  
- **--length**: Sequence length bounds. Provide one value (upper bound) or two values (min and max).  
- **--quality**: Minimum average quality threshold.

**Example usage:**

```bash
python fastq_filter_tool.py --input example.fastq --output filtered.fastq --gc 40 60 --length 50 150 --quality 30
```
---
## Logging
The tool now also supports logging using Python's built-in logging module. All log messages are written to the file fastq_filter.log.

**Informational messages**: The tool logs an informational message on start-up and after successful filtering.

**Error logging**: If an exception occurs during filtering, an error message is logged.

---
## Testing
A total of 8 tests cover the following aspects:

**Error Handling**: At least one test verifies that invalid GC bounds or a non-existent input file results in an error.

**File I/O Operations**: Tests ensure that the output file is created and correctly written, and that the log file is generated.

**CLI Argument Parsing**: Tests confirm that the CLI tool returns an error when required arguments are missing.

The tests are organized in test_fastq_filter.py and can be run from the root directory of the project using:
```bash
python -m unittest discover
```

## Contact
For any questions or feedback please contact:

Elena Dementeva
Email: dementevayelena@gmail.com

