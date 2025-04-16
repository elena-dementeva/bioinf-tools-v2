from typing import Tuple, Union
from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


class BiologicalSequence(ABC):
    """
    Abstract class for biological sequences.
    """

    def __init__(self, sequence: str):
        self.sequence = sequence.upper()
        self.check_alphabet()

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.sequence}')"

    @abstractmethod
    def check_alphabet(self):
        """ Validates the sequence alphabet. """
        pass

    def to_oneline_fasta(self, output_fasta: str) -> None:
        """
        Writes the sequence in single-line FASTA format.
        """
        with open(output_fasta, 'w') as outfile:
            outfile.write(f">{self.__class__.__name__}\n{self.sequence}\n")


class NucleicAcidSequence(BiologicalSequence):
    """
    Base class for nucleic acid sequences (DNA and RNA).
    """

    complement_map = {}

    def __init__(self, sequence: str):
        if not self.complement_map:
            raise NotImplementedError("NucleicAcidSequence is an abstract class and cannot be instantiated directly.")
        super().__init__(sequence)

    def check_alphabet(self):
        """ Validates that the sequence contains only valid nucleotides. """
        valid_nucleotides = set(self.complement_map.keys())
        if not set(self.sequence).issubset(valid_nucleotides):
            raise ValueError(f"Invalid sequence: {self.sequence}")

    def complement(self):
        """ Returns the complement of the sequence. """
        return self.__class__(''.join(self.complement_map[nuc] for nuc in self.sequence))

    def reverse(self):
        """ Returns the reversed sequence. """
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self):
        """ Returns the reverse complement of the sequence. """
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    """
    Class for DNA sequences.
    """

    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    def transcribe(self):
        """ Converts DNA to RNA. """
        return RNASequence(self.sequence.replace('T', 'U'))


class RNASequence(NucleicAcidSequence):
    """
    Class for RNA sequences.
    """

    complement_map = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}


class AminoAcidSequence(BiologicalSequence):
    """
    Class for amino acid sequences (proteins).
    """

    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

    def check_alphabet(self):
        """ Validates that the sequence contains only valid amino acids. """
        if any(aa not in self.valid_amino_acids for aa in self.sequence):
            raise ValueError(f"Invalid amino acid sequence: {self.sequence}")

    def count_hydrophobic_residues(self):
        """ Counts the number of hydrophobic residues in the sequence. """
        hydrophobic = set("AVLIMFWP")
        return sum(1 for aa in self.sequence if aa in hydrophobic)


class BlastParser:
    """
    Class for parsing BLAST output files.
    """

    @staticmethod
    def parse_blast(input_file: str, output_file: str) -> None:
        """
        Extracts the best match descriptions from a BLAST output file and saves them alphabetically.
        """
        descriptions = set()
        with open(input_file) as infile:
            for line in infile:
                if "Sequences producing significant alignments:" in line:
                    break
            for line in infile:
                if not line.strip():
                    break
                parts = line.split(None, 1)
                if len(parts) > 1:
                    descriptions.add(parts[1])

        with open(output_file, 'w') as outfile:
            for desc in sorted(descriptions):
                outfile.write(f"{desc}\n")


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[Tuple[float, float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0
) -> None:
    """
    Reads a FASTQ file, filters sequences based on GC content, length, quality,
    and writes the filtered sequences to a new file.
    """

    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)

    filtered_records = []

    for record in SeqIO.parse(input_fastq, "fastq"):
        gc_content = gc_fraction(record.seq) * 100
        seq_len = len(record)
        avg_quality = sum(record.letter_annotations["phred_quality"]) / len(record)

        if (
            (gc_bounds[0] <= gc_content <= gc_bounds[1])
            and (length_bounds[0] <= seq_len <= length_bounds[1])
            and avg_quality >= quality_threshold
        ):
            filtered_records.append(record)

    SeqIO.write(filtered_records, output_fastq, "fastq")
