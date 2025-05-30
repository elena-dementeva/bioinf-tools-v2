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
cd bioinf-tools
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

## Contact
For any questions or feedback please contact:

Elena Dementeva
Email: dementevayelena@gmail.com

