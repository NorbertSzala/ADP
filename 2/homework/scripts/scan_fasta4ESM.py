#!/usr/bin/python3

"""
scan_fasta4ESM.py

Scan a directory for FASTA files that have not yet been processed by predictESM.py.

A FASTA file is considered processed if a directory with the same name
(without .fas extension) already exists.

Example:
fastas/Mycgen/Mycgen_000-019.fas  ->  pdb/Mycgen/Mycgen_000-019/

If the directory does not exist, the script runs predictESM.py for that file.

Usage:
python3 scan_fasta4ESM.py ./fastas/

Important:
Only ONE FASTA file is processed per execution.
"""

import argparse
from pathlib import Path
import subprocess
import sys


################################
# Argument parsing
################################


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="scan_fasta4ESM.py",
        description="Scan folder for FASTA files not yet processed by ESM",
        usage="%(prog)s -I/--input [options]",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-I", "--input", help="Directory containing FASTA files", type=Path
    )

    args = parser.parse_args()

    if not args.input.exists():
        parser.error("Provided path does not exist")

    if not args.input.is_dir():
        parser.error("Provided path is not a directory")

    return args


################################
# Main
################################


def main():

    args = parse_arguments()
    fasta_dir = args.input

    fasta_files = sorted(fasta_dir.rglob("*.fas"))

    for fasta in fasta_files:
        dataset = fasta.parent.name
        range_part = fasta.stem.split("_")[-1]
        output_dir = Path(f"results/pdb/{dataset}/{dataset}_{range_part}")
        if output_dir.exists():
            continue

        print(f"Processing {fasta}")
        output_prefix = f"results/pdb/{dataset}/{dataset}_"

        subprocess.run(
            [
                "python3",
                "scripts/predictESM.py",
                "-I",
                str(fasta),
                "-O",
                output_prefix,
            ]
        )

        break


if __name__ == "__main__":
    main()
