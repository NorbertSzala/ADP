#!/usr/bin/python3

"""
fasta_splitter.py

Description
-----------
Split a FASTA file into multiple smaller FASTA files.

The script:
- reads a FASTA file containing many protein sequences
- filters sequences longer than a specified maximum length
- splits remaining sequences into files containing a fixed number of sequences

Each output file contains at most S sequences.

Usage
-----
python3 fasta_splitter.py -I input -O ouptut -M max_length -S max_sequences

Arguments
---------
-I / --input
    Input FASTA file (.fa, .fasta, .fas)

-O / --output
    Output prefix including directory and file prefix.
    Example: ./results/Mycgen_

-M / --max_length
    Maximum allowed sequence length (default: 400)

-S / --max_sequences
    Maximum number of sequences per output FASTA file (default: 20)

Example
-------
python3 ./scripts/fasta_splitter.py \
-I ./data/Mycoplasmoides_genitalium_G37.fas \
-O ./results/fastas/Mycgen_ \
-M 400 \
-S 20

Output
------
./results/fastas/Mycgen_000-019.fas
./results/fastas/Mycgen_020-039.fas
...
./results/fastas/Mycgen_420-433.fas
"""

import argparse
from pathlib import Path
import sys


################################
# Helper functions
################################


def parse_arguments():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser(
        prog="fasta_splitter.py",
        description="Script that takes a FASTA file, an output_name, max_size, and size as input to produce the proper files",
        usage="%(prog)s -I/--input [options]",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-I",
        "--input",
        help="Input fasta file (allowed different extensions, .fa, .fasta, .fas)",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "-O",
        "--output",
        help="Path to folder AND PREFFIX where to save splitted fasta files. In example: ./results/prefix_",
        required=False,
        type=Path,
        default=None,
    )

    parser.add_argument(
        "-M",
        "--max_length",
        help="Maximal length of parsed protein sequences from the input file",
        type=int,
        default=400,
    )

    parser.add_argument(
        "-S",
        "--max_sequences",
        help="Maximal number of sequences in a single .fasta file",
        type=int,
        default=20,
    )
    return parser.parse_args()


def create_output_file(input_file: Path, output_file: Path = None) -> Path:
    if output_file is None:
        output_dir = input_file.parent / "output"
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / f"{input_file.stem}.fas"
    else:
        output_file.parent.mkdir(parents=True, exist_ok=True)
    return output_file


def parse_fasta_records(test: str):
    """yield fasta records as (header, sequence)"""
    header, seq_lines = None, []

    for line in test.splitlines():
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


################################
# Main functions
################################


def main():
    ######### read arguments #########
    print(
        f'Running {Path(__file__).name} with arguments: {" ".join(map(str, sys.argv[1:]))}'
    )
    args = parse_arguments()
    input_file = args.input
    if not input_file.exists() or input_file.suffix not in (".fas", ".fasta", ".fa"):
        raise ValueError(
            "Input file should be saved in .fasta format (also acceptable .fa, .fas) or does not exists in given location"
        )

    output_file = create_output_file(input_file, args.output)
    max_len = args.max_length
    max_seq = args.max_sequences
    if max_seq <= 0:
        raise ValueError("max_sequences must be greater than 0")
    if max_len <= 0:
        raise ValueError("max_length must be greater than 0")
    records = []

    ######### Read input file and create list with sequences #########
    with open(input_file) as f:
        for header, seq in parse_fasta_records(f.read()):

            seq = seq.strip()

            if len(seq) <= max_len:  # reject proteins longer than given argument
                records.append((header, seq))
    output_prefix = str(output_file)

    ######### Save results #########
    for start in range(0, len(records), max_seq):

        end = min(start + max_seq, len(records))

        chunk = records[start:end]

        file_name = f"{output_prefix}{start:03d}-{end-1:03d}.fas"
        file_path = Path(file_name)

        file_path.parent.mkdir(parents=True, exist_ok=True)

        with open(file_path, "w") as out:
            for header, seq in chunk:
                out.write(f"{header}\n{seq}\n")


if __name__ == "__main__":
    main()
