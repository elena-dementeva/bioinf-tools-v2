from typing import Dict, Tuple, Union
from scripts.dna_rna_tools import (
    transcribe,
    reverse,
    complement,
    reverse_complement,
    validate_sequence
)

from scripts.fastq_utils import read_fastq, write_fastq

def run_dna_rna_tools(*args: str) -> Union[str, list]:
    """
    Processes DNA/RNA sequences based on the specified action.

    Args:
        *args: Sequences followed by an action ("transcribe", "reverse",
        "complement" or "reverse_complement").

    Returns:
        Processed sequences as a string or a list of strings.
    """
    *seqs, action = args

    results = []

    for seq in seqs:
        validate_sequence(seq)
        if action == "transcribe":
            results.append(transcribe(seq))
        elif action == "reverse":
            results.append(reverse(seq))
        elif action == "complement":
            results.append(complement(seq))
        elif action == "reverse_complement":
            results.append(reverse_complement(seq))
        else:
            raise ValueError("error")

    return results[0] if len(results) == 1 else results


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[Tuple[float, float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0
) -> None:
    """
    Reads a FASTQ file, filters sequences based on GC content, length, and quality, 
    and writes the filtered sequences to a new FASTQ file.

    Args:
        input_fastq: Path to the input FASTQ file.
        output_fastq: Path for saving filtered sequences.
        gc_bounds: Tuple or value which specify GC content bounds.
        length_bounds: Tuple or value which specify length bounds.
        quality_threshold: Minimum quality for filtering.
    """
    with open(output_fastq, 'w') as out_file:
        for name, seq, qual in read_fastq(input_fastq):
            gc_content = (sum(base in 'GCgc' for base in seq) / len(seq)) * 100
            seq_len = len(seq)
            avg_quality = sum(ord(char) - 33 for char in qual) / len(qual)

            if (
                (gc_bounds[0] <= gc_content <= gc_bounds[1]) and
                (length_bounds[0] <= seq_len <= length_bounds[1]) and
                avg_quality >= quality_threshold
            ):
                write_fastq(out_file, name, seq, qual)
