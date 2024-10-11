def read_fastq(input_fastq):
    """
    Read sequences from a FASTQ file.

    Args:
        input_fastq: Path to FASTQ file.
    """
    with open(input_fastq, 'r') as file:
        while True:
            header = file.readline().strip()
            if not header:
                break
            sequence = file.readline().strip()
            file.readline()  # Skip '+'
            quality = file.readline().strip()
            yield header, sequence, quality

def write_fastq(output_fastq, header, sequence, quality):
    """
    Writes a sequence to a FASTQ file.

    Args:
        output_fastq (str): Output FASTQ file path.
        header (str): Sequence header.
        sequence (str): DNA sequence.
        quality (str): Quality scores.
    """
    with open(output_fastq, 'a') as file:
        file.write(f"{header}\n{sequence}\n+\n{quality}\n")
