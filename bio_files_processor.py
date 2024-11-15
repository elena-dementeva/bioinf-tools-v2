def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str = "output_oneline.fasta"
) -> None:
    """
    Converts a multiline FASTA file to single-line format for each sequence.
    """
    with open(input_fasta) as infile, open(output_fasta, 'w') as outfile:
        sequence, header = "", None
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    outfile.write(f"{header}\n{sequence}\n")
                header, sequence = line, ""
            else:
                sequence += line
        if header:
            outfile.write(f"{header}\n{sequence}\n")


def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    Extracts best match descriptions from a BLAST output file
        and saves them alphabetically.
    """
    descriptions = set()
    with open(input_file) as infile:
        for line in infile:
            if "Sequences producing significant alignments:" in line:
                break
        for line in infile:
            if not line.strip():
                break
            description = line.split(None, 1)[1]
            if description:
                descriptions.add(description)
    with open(output_file, 'w') as outfile:
        for desc in sorted(descriptions):
            outfile.write(f"{desc}\n")
