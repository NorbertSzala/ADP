#!/usr/bin/python3

"""
Description:
gp2fasta.py converts GenPept (.gp) protein files into FASTA (.fas) format.

The script:
- works with multiple records in one file
- extracts protein sequences from the ORIGIN section
- builds FASTA headers
- can modify organism name format
- can filter sequences by length
- can include Locus or GI in the header

Usage:
./scripts/gp2fasta.py -I ./data/input.gp

Additional options:
Convert organism name format:
-C 1   # Homo sapiens
-C 2   # H. sapiens
-C 3   # Homsap

Choose ID type:
--id_type Locus
--id_type GI

Add extra fields to FASTA header:
-L      # add locus
-G      # add GI

Filter by length:
--length_min 100
--length_max 300

Example
./scripts/gp2fasta.py -I ./data/example.gp -C 2 --id_type GI -G --length_min 100 -O results/output.fas
"""
import argparse
from pathlib import Path
import re
import sys


################################
# Helper functions
################################


def parse_arguments():
    """Parse command-line arguments for gp2fasta.py"""
    parser = argparse.ArgumentParser(
        prog="gp2fasta.py",
        description="Program converting GeneBank format (.gp) to Fasta format (.fas)",
        usage="%(prog)s -I/--input [options]",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-I", "--input", help="Input GeneBank file (.gp)", required=True, type=Path
    )
    parser.add_argument(
        "-O",
        "--output",
        help="Output Fasta file (.fas)",
        required=False,
        type=Path,
        default=None,
    )
    parser.add_argument(
        "-C",
        "--convert_names",
        help="Convert names to standard format, give option 1, 2 or 3. When 1 is given, it converts the  (Homo Sapiens) to Homo sapiens. When 2 is given, it converts the name to H. sapiens. When 3 is given, it converts the name to Homsap.\n It is not possible to convert the name from shorter to longer format",
        required=False,
        default=1,
        choices=[1, 2, 3],
        type=int,
    )
    parser.add_argument(
        "-L",
        "--extract_locus",
        help="Extract locus name from GeneBank file",
        required=False,
        action="store_true",
    )
    parser.add_argument(
        "-G",
        "--extract_GI",
        help="Extract GI name from GeneBank file",
        required=False,
        action="store_true",
    )
    parser.add_argument(
        "--id_type",
        help="Extract ID from GeneBank file. Specify whether you want to have Locus (default) or GI",
        required=False,
        choices=["Locus", "GI"],
        default="Locus",
        type=str,
    )
    parser.add_argument(
        "--length_min",
        help="Minimum sequence length",
        required=False,
        type=int,
    )
    parser.add_argument(
        "--length_max",
        help="Maximum sequence length",
        required=False,
        type=int,
    )
    parser.add_argument(
        "-D",
        "--delimiter",
        help="Delimiter used in FASTA header",
        default="|",
        type=str,
    )
    parser.add_argument(
        "-R",
        "--rawheader",
        help="Leave only protein ID in FASTA header (ignore other header options)",
        action="store_true",
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


def read_gp_file(gp_file):
    with open(gp_file, "r") as f:
        content = f.read()
    return content


def parse_gp_records(text: str):
    """Split the input file into records by the separator // and return each record as text."""
    for block in text.split("\n//"):
        if block.strip():
            yield block.strip()


def parse_record(block: str, idtype: str, name_option: int):
    """
    Parse a single record and return a dict:
    {
      "id": ...,
      "organism": ...,
      "length": int,
      "seq": "...."
    }
    """
    locus_len = None
    version_line = None
    organism = None
    seq_lines = []
    in_origin = False

    for line in block.splitlines():
        if line.startswith("LOCUS"):
            parts = line.split()
            # LOCUS <locus> <length> aa ...
            if len(parts) >= 3 and parts[2].isdigit():
                locus_len = int(parts[2])

        elif line.startswith("VERSION"):
            version_line = line

        elif line.startswith("  ORGANISM"):
            # "   ORGANISM Homo sapiens"
            organism = line.replace("  ORGANISM ", "").strip()

        elif line.startswith("ORIGIN"):
            in_origin = True

        elif in_origin:
            # Line with sequence data, take only protein sequence
            letters = re.sub(r"[^A-Za-z]", "", line)
            if letters:
                seq_lines.append(letters.upper())
    seq = "".join(seq_lines)

    # id from version (locus or GI depending on --id_type)
    rec_id = id_type(version_line, idtype) if version_line else None

    # always try to capture GI (optional in header)
    gi = extract_gi_from_version_line(version_line) if version_line else None

    org_conv = convert_names(organism, name_option) if organism else None

    return {
        "id": rec_id,
        "gi": gi,
        "organism": org_conv,
        "length": len(seq),
        "seq": seq,
    }


def passes_length_filter(
    length: int | None, seq: str, minlen: int | None, maxlen: int | None
) -> bool:
    """
    If length is given in header, use it. Otherwise, calculate from sequence.
    """
    if length is None:
        length = len(seq)

    if minlen is not None and length < minlen:
        return False
    if maxlen is not None and length > maxlen:
        return False
    return True


def record_to_fasta(
    rec: dict,
    extract_locus_flag: bool,
    extract_gi_flag: bool,
    delimiter: str,
    rawheader=bool,
) -> str:
    """
    Convert a record dict to FASTA format string.
    The header includes ID, organism, length, and optionally locus name if extract_locus_flag is True.
    """
    rec_id = rec["id"] or "UNKNOWN_ID"

    # RAW HEADER MODE
    if rawheader:
        header = f">{rec_id}"
    else:
        org = rec["organism"] or "UNKNOWN_ORGANISM"
        length = rec["length"] if rec["length"] is not None else len(rec["seq"])

        header_parts = [rec_id, org, f"len={length}"]

        if extract_locus_flag:
            header_parts.append(f"locus={rec_id}")

        if extract_gi_flag and rec.get("gi"):
            header_parts.append(f"GI:{rec['gi']}")

        header = ">" + delimiter.join(header_parts)
    return header + "\n" + str(rec["seq"]).upper() + "\n"


def convert_names(name: str, option: int) -> str:
    """
    Convert organism name to selected format.

    Option 1: Homo sapiens
    Option 2: H. sapiens
    Option 3: Homsap

    If format is not 'Genus species', return original.
    """

    name = name.strip()

    match = re.fullmatch(r"([A-Z][a-z]+) ([a-z]+)", name)

    if not match:
        return name

    genus, species = match.groups()

    if option == 1:
        return f"{genus} {species}"
    elif option == 2:
        return f"{genus[0]}. {species}"
    elif option == 3:
        return f"{genus[:3]}{species[:3]}"

    return name


def id_type(line: str, locus_type: str) -> str:
    """
    Extract ID from VERSION line.

    Example:
    VERSION     XP_054228232.1 GI:123234345

    returns XP_054228232 if locus_type == Locus, returns 123234345 if locus_type == GI
    """

    parts = line.strip().split()

    if locus_type == "Locus":
        # XP_054228232.1 - remove version (.1)
        return parts[1].split(".")[0]

    elif locus_type == "GI":
        match = re.search(r"GI:(\d+)", line)
        if match:
            return match.group(1)

    return parts[1].split(".")[0]


def extract_gi_from_version_line(line: str) -> str | None:
    """Extract GI number from VERSION line, e.g. '... GI:123' -> '123'."""
    if not line:
        return None
    m = re.search(r"GI:(\d+)", line)
    return m.group(1) if m else None


################################
# Main functions
################################


def main():
    # read arguments
    print(
        f'Running {Path(__file__).name} with arguments: {" ".join(map(str, sys.argv[1:]))}'
    )
    # parse arguments, validate input file, create output file path, and read input file
    args = parse_arguments()
    input_file = args.input
    if not input_file.exists() or input_file.suffix != ".gp":
        raise ValueError("Input file must exist and have .gp extension")
    output_file = create_output_file(input_file, args.output)
    convert_names_option = args.convert_names
    extract_locus_option = args.extract_locus
    extract_gi = args.extract_GI
    idtype = args.id_type
    lmin = args.length_min
    lmax = args.length_max
    delimiter = args.delimiter
    rawheader = args.rawheader
    if lmin is not None and lmax is not None and lmin > lmax:
        raise ValueError("Minimum length cannot be greater than maximum length")

    # Main processing: read input file, parse records, filter by length, and save to output file
    records = []

    with open(input_file) as f:
        for block in parse_gp_records(f.read()):
            # read input file   # iterate over records and extract data when option in argauments is given
            rec = parse_record(block, idtype, convert_names_option)

            # if sequence is empty, skip record
            if not rec["seq"]:
                continue

            if passes_length_filter(rec["length"], rec["seq"], lmin, lmax):
                records.append(rec)

    # save records to output file in FASTA format
    with open(output_file, "w") as out:
        for rec in records:
            out.write(
                record_to_fasta(
                    rec, extract_locus_option, extract_gi, delimiter, rawheader
                )
            )

    print(f"Conversion completed. Output saved to {output_file}")


if __name__ == "__main__":
    main()
