def transcribe(sequence: str) -> str:
    # РНК не должна быть для транскрипции
    if 'U' in sequence or 'u' in sequence:
        raise ValueError("RNA detected. Only DNA can be transcribed.")
    return sequence.replace('T', 'U').replace('t', 'u')


def reverse(sequence: str) -> str:
    return sequence[::-1]


def complement(sequence: str) -> str:
    complement_dict = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
        'U': 'A', 'u': 'a'
    }
    return ''.join(complement_dict[nuc] for nuc in sequence)


def reverse_complement(sequence: str) -> str:
    return reverse(complement(sequence))


def validate_sequence(sequence: str) -> None:
    valid_dna = "ATCGatcg"
    valid_rna = "AUCGaucg"
    # T и U не должны быть вместе
    if 'U' in sequence or 'u' in sequence:
        if 'T' in sequence or 't' in sequence:
            raise ValueError("contains both T and U")
        if any(nucleotide not in valid_rna for nucleotide in sequence):
            raise ValueError("Invalid for RNA")
    else:
        if any(nucleotide not in valid_dna for nucleotide in sequence):
            raise ValueError("Invalid for DNA")


def run_dna_rna_tools(*args: str):
    # Распаковка аргументов
    *sequences, procedure = args

    results = []

    for seq in sequences:
        validate_sequence(seq)
        if procedure == "transcribe":
            results.append(transcribe(seq))
        elif procedure == "reverse":
            results.append(reverse(seq))
        elif procedure == "complement":
            results.append(complement(seq))
        elif procedure == "reverse_complement":
            results.append(reverse_complement(seq))
        else:
            raise ValueError("Procedure is unknown")

    # если одна последовательность -строкf, несколько - список
    return results[0] if len(results) == 1 else results
