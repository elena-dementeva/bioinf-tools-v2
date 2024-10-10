from typing import Dict, Tuple, Union
from scripts.dna_rna_tools import (
    transcribe,
    reverse,
    complement,
    reverse_complement,
    validate_sequence
)


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
    seqs: Dict[str, Tuple[str, str]],
    gc_bounds: Union[Tuple[float, float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0
) -> Dict[str, Tuple[str, str]]:
    """    
    Filters sequences based on GC content, length and quality.

    Args:
        seqs: dictionary with sequence names as keys and 
        quality as value.
        gc_bounds: Tuple or value which specify GC content bounds.
        length_bounds: Tuple or value which specify length bounds.
        quality_threshold: Minimum quality for filtering.

    Returns:
        Filtered sequences as a dictionary.
    """

    def calc_gc(seq: str) -> float:
        gc_count = sum(base in 'GCgc' for base in seq)
        return (gc_count / len(seq)) * 100 if seq else 0

    def calc_quality(qual: str) -> float:
        return sum(ord(char) - 33 for char in qual) / len(qual) if qual else 0

    def check_bounds(
        val: Union[int, float],
        bounds: Union[
            Tuple[Union[int, float], Union[int, float]], Union[int, float]
        ]
    ) -> bool:
        if isinstance(bounds, tuple):
            return bounds[0] <= val <= bounds[1]
        return val <= bounds

    if isinstance(gc_bounds, (float, int)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, (int, float)):
        length_bounds = (0, length_bounds)

    filtered = {}

    for name, (seq, qual) in seqs.items():
        gc_content = calc_gc(seq)
        seq_len = len(seq)
        avg_quality = calc_quality(qual)

        if (
            check_bounds(gc_content, gc_bounds) and
            check_bounds(seq_len, length_bounds) and
            avg_quality >= quality_threshold
        ):
            filtered[name] = (seq, qual)

    return filtered
