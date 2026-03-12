#!/usr/bin/python3

"""
predictESM.py

Description
-----------
Simple wrapper that sends protein sequences from a multi-sequence FASTA file
to the ESMFold API and saves predicted structures.

Usage
-----
python3 predictESM.py ./fastas/Mycgen_000-019.fas

Output
------
Creates a directory named after the FASTA file and stores predicted structures:

./fastas/Mycgen_000-019/WP_009885972_esmfold_v1.pdb
./fastas/Mycgen_000-019/WP_009885652_esmfold_v1.pdb
./fastas/Mycgen_000-019/WP_009885673_esmfold_v1.pdb

Notes
-----
• Sequences are sent one by one to the ESMFold API
• 5 second delay between requests
• Maximum sequence length allowed by API is ~400 aa
"""

import argparse
from pathlib import Path
import sys
import time
import requests


################################
# Helper functions
################################
def parse_arguments():
    """Parse command line arguments for predictESM.py"""

    parser = argparse.ArgumentParser(
        prog="predictESM.py",
        description="Send protein sequences from a FASTA file to the ESMFold API and store predicted structures.",
        usage="%(prog)s -I/--input [options]",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-I",
        "--input",
        help="Input FASTA file",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "-O",
        "--output",
        help="Output prefix directory (example: results/pdb/Mycgen_)",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "-d",
        "--delay",
        help="Delay between API calls in seconds",
        type=int,
        default=5,
    )

    args = parser.parse_args()

    # validate input
    if not args.input.exists():
        parser.error("Input FASTA file does not exist")

    if args.input.suffix not in (".fas", ".fasta", ".fa"):
        parser.error("Input file must be FASTA (.fas, .fasta, .fa)")

    if args.delay < 0:
        parser.error("Delay must be >= 0")

    return args


def parse_fasta_records(text):
    """Yield FASTA records as (header, sequence)"""

    header = None
    seq_lines = []

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header:
                yield header, "".join(seq_lines)

            header = line
            seq_lines = []

        else:
            seq_lines.append(line)

    if header:
        yield header, "".join(seq_lines)


def get_protein_id(header):
    """Extract protein ID from FASTA header"""
    return header.split()[0].replace(">", "")


def send_to_esm(sequence, retries=3, wait=10):
    """Send sequence to ESMFold API with retry logic"""
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    for attempt in range(retries):
        try:
            response = requests.post(url, data=sequence, timeout=120)
            if response.status_code == 200:
                return response.text
            print(f"API error {response.status_code}, retry {attempt+1}/{retries}")
        except requests.exceptions.RequestException as e:
            print(f"Request failed: {e}")
        time.sleep(wait)

    print("Failed after retries")
    return None


################################
# Main functions
################################


def main():

    print(
        f'Running {Path(__file__).name} with arguments: {" ".join(map(str, sys.argv[1:]))}'
    )

    args = parse_arguments()
    input_file = args.input
    output_prefix = args.output
    delay = args.delay

    # extract numeric range from filename
    # Mycgen_000-019.fas -> 000-019
    range_part = input_file.stem.split("_")[-1]
    output_dir = Path(f"{output_prefix}{range_part}")

    output_dir.mkdir(parents=True, exist_ok=True)

    with open(input_file) as f:
        for header, seq in parse_fasta_records(f.read()):
            seq = seq.strip()

            if len(seq) > 400:
                print(f"Skipping {header} (sequence longer than 400 aa)")
                continue

            protein_id = get_protein_id(header)

            output_file = output_dir / f"{protein_id}_esmfold_v1.pdb"

            if output_file.exists():
                print(f"Skipping {protein_id} (already predicted)")
                continue

            pdb_structure = send_to_esm(seq)

            if pdb_structure is None:
                print(f"Skipping {protein_id} after repeated failures")
                continue

            output_file = output_dir / f"{protein_id}_esmfold_v1.pdb"
            print(f"Predicting {protein_id} -> {output_file}")

            with open(output_file, "w") as out:
                out.write(pdb_structure)

            time.sleep(delay)


if __name__ == "__main__":
    main()
