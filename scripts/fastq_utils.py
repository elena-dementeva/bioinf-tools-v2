import os


def read_fastq(file_path: str) -> dict:
    """
    Reads sequences from a FASTQ file and returns them as a dictionary.
    Key is the sequence name, value is a tuple (sequence, quality).
    """
    sequences = {}
    with open(file_path, 'r') as file:
        while True:
            name = file.readline().strip()
            if not name:
                break
            seq = file.readline().strip()
            file.readline()
            qual = file.readline().strip()
            sequences[name] = (seq, qual)
    return sequences


def write_fastq(sequences: dict, output_fastq: str) -> None:
    """
    Writes a sequence to a FASTQ file.
    """
    if not output_fastq.endswith(".fastq"):
        output_fastq += ".fastq"
    output_path = os.path.join("filtered", output_fastq)
    output_dir = os.path.dirname(output_path)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if os.path.exists(output_path):
        raise FileExistsError("File exists. Choose a different name.")

    with open(output_path, 'w') as file:
        for name, (seq, qual) in sequences.items():
            file.write("{}\n{}\n+\n{}\n".format(name, seq, qual))
