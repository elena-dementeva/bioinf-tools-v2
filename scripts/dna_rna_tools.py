def transcribe(sequence: str) -> str:
    """
    Converts a DNA sequence to RNA by replacing 'T' with 'U'.
    Args:
        sequence (str): DNA sequence.
    Returns:
        str: RNA sequence.
    Raises:
        ValueError: If the sequence contains RNA ('U' or 'u').
    """
    if 'U' in sequence or 'u' in sequence:
        raise ValueError("RNA detected. Only DNA can be transcribed.")
    return sequence.replace('T', 'U').replace('t', 'u')


def reverse(sequence: str) -> str:
    """
    Reverse the sequence.
    Args:
        sequence (str): Input sequence.
    Returns:
        str: Reversed sequence.
    """
    return sequence[::-1]


def complement(sequence: str) -> str:
    """
    Returns the complement of a DNA/RNA sequence.
    Args:
        sequence (str): Input sequence.
    Returns:
        str: Complement sequence.
    """
    complement_dict = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
        'U': 'A', 'u': 'a'
    }
    return ''.join(complement_dict[nuc] for nuc in sequence)


def reverse_complement(sequence: str) -> str:
    """
    Returns the reverse complement of a sequence.
    Args:
        sequence (str): Input sequence.
    Returns:
        str: Reverse complement sequence.
    """
    return reverse(complement(sequence))


def validate(sequence: str) -> None:
    """
    Checks if a sequence is DNA or RNA.
    Args:
        sequence (str): Sequence to validation.
    Raises:
        ValueError: If the sequence is invalid or contains T and U.
    """
    valid_dna = "ATCGatcg"
    valid_rna = "AUCGaucg"
    if 'U' in sequence or 'u' in sequence:
        if 'T' in sequence or 't' in sequence:
            raise ValueError("contains both T and U")
        if any(nucleotide not in valid_rna for nucleotide in sequence):
            raise ValueError("Invalid for RNA")
    else:
        if any(nucleotide not in valid_dna for nucleotide in sequence):
            raise ValueError("Invalid for DNA")
