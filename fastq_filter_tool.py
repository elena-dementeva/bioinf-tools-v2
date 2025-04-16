import argparse
import logging
from bioinf_tools import filter_fastq

def parse_args():
    """
    Parses command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="FastQ Filter Tool")
    parser.add_argument("--input", required=True, help="Path to the input FASTQ file")
    parser.add_argument("--output", required=True, help="Path to the output FASTQ file")
    parser.add_argument("--gc", type=float, nargs='+', default=[0, 100],
                        help="GC content bounds: specify one value (upper bound) or two values (min and max)")
    parser.add_argument("--length", type=int, nargs='+', default=[0, 1000000],
                        help="Sequence length bounds: specify one value (upper bound) or two values (min and max)")
    parser.add_argument("--quality", type=float, default=0,
                        help="Minimum average quality threshold for sequences")
    return parser.parse_args()

def setup_logging():
    """
    Configures logging to output messages to a file.
    """
    logging.basicConfig(
        filename="fastq_filter.log",
        level=logging.DEBUG,
        format="%(asctime)s %(levelname)s: %(message)s"
    )

def main():
    """
    Main entry point for the FastQ Filter Tool.
    
    It parses command-line arguments, sets up logging, and calls the filtering function.
    Logs an informational message when starting and an error message if an exception occurs.
    """
    args = parse_args()
    setup_logging()
    logging.info("Starting FastQ Filter Tool with arguments: %s", args)

    try:
        # Process GC bounds arguments: if one value provided, assume lower bound is 0.
        if len(args.gc) == 1:
            gc_bounds = (0, args.gc[0])
        elif len(args.gc) == 2:
            gc_bounds = tuple(args.gc)
        else:
            raise ValueError("Invalid number of GC bounds provided.")
        
        # Process length bounds arguments: if one value provided, assume lower bound is 0.
        if len(args.length) == 1:
            length_bounds = (0, args.length[0])
        elif len(args.length) == 2:
            length_bounds = tuple(args.length)
        else:
            raise ValueError("Invalid number of length bounds provided.")

        filter_fastq(args.input, args.output, gc_bounds, length_bounds, args.quality)
        logging.info("FASTQ filtering completed successfully.")
    except Exception as e:
        logging.error("An error occurred during FASTQ filtering: %s", e)
        raise

if __name__ == "__main__":
    main()
