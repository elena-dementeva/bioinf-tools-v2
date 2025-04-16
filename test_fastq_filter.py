import os
import tempfile
import unittest
import subprocess
from bioinf_tools import filter_fastq

class TestFastqFilterTool(unittest.TestCase):
    """Test suite for FastQ filtering functionality and CLI behavior."""

    def setUp(self):
        """
        Sets up temporary input and output FASTQ files for testing.
        Creates a temporary input FASTQ file with two example records.
        """
        self.input_fastq = tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fastq")
        self.input_fastq.write("@seq1\nACGT\n+\n!!!!\n")
        self.input_fastq.write("@seq2\nAAAA\n+\n!!!!\n")
        self.input_fastq.close()
        self.output_fastq = tempfile.NamedTemporaryFile(delete=False, suffix=".fastq").name

    def tearDown(self):
        """
        Cleans up temporary files created during testing.
        """
        if os.path.exists(self.input_fastq.name):
            os.remove(self.input_fastq.name)
        if os.path.exists(self.output_fastq):
            os.remove(self.output_fastq)
        if os.path.exists("fastq_filter.log"):
            os.remove("fastq_filter.log")

    def test_output_file_creation(self):
        """
        Test that the filter_fastq function creates an output file.
        """
        filter_fastq(self.input_fastq.name, self.output_fastq)
        self.assertTrue(os.path.exists(self.output_fastq), "Output file was not created.")

    def test_output_file_not_empty(self):
        """
        Test that the output file is not empty when sequences pass the filters.
        """
        filter_fastq(self.input_fastq.name, self.output_fastq, gc_bounds=(0, 100), length_bounds=(0, 10), quality_threshold=0)
        with open(self.output_fastq, "r") as f:
            content = f.read()
        self.assertNotEqual(content, "", "Output file is empty despite expected sequences.")

    def test_invalid_gc_bounds_error(self):
        """
        Test that an error is raised when GC bounds are invalid (min > max).
        """
        with self.assertRaises(Exception):
            filter_fastq(self.input_fastq.name, self.output_fastq, gc_bounds=(50, 40))

    def test_length_filtering(self):
        """
        Test that filtering by length works as expected: if upper bound is too small, no sequences pass.
        """
        filter_fastq(self.input_fastq.name, self.output_fastq, length_bounds=(0, 3))
        with open(self.output_fastq, "r") as f:
            content = f.read()
        self.assertEqual(content, "", "Output file should be empty due to strict length filtering.")

    def test_quality_filtering(self):
        """
        Test that a high quality threshold filters out records.
        """
        # Create a temporary FASTQ file with a single record.
        with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fastq") as tmp:
            tmp.write("@seq3\nACGT\n+\n!!!!\n")
            tmp_name = tmp.name
        try:
            filter_fastq(tmp_name, self.output_fastq, quality_threshold=40)
            with open(self.output_fastq, "r") as f:
                content = f.read()
            self.assertEqual(content, "", "Output file should be empty when quality threshold is high.")
        finally:
            os.remove(tmp_name)

    def test_exception_for_nonexistent_input(self):
        """
        Test that filter_fastq raises an exception when the input file does not exist.
        """
        with self.assertRaises(Exception):
            filter_fastq("nonexistent.fastq", self.output_fastq)

    def test_logging_file_creation(self):
        """
        Test that running the CLI script creates a log file.
        """
        cmd = ["python", "fastq_filter_tool.py", "--input", self.input_fastq.name, "--output", self.output_fastq]
        subprocess.call(cmd)
        self.assertTrue(os.path.exists("fastq_filter.log"), "Log file was not created by the CLI tool.")

    def test_argparse_error_on_missing_required_args(self):
        """
        Test that the CLI script fails when required command-line arguments are missing.
        """
        cmd = ["python", "fastq_filter_tool.py", "--input", self.input_fastq.name]
        result = subprocess.run(cmd, capture_output=True, text=True)
        self.assertIn("error", result.stderr.lower(), "CLI tool did not report error when required arguments were missing.")

if __name__ == "__main__":
    unittest.main()
