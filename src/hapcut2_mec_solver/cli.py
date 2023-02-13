import argparse
import sys
import os
import json
from typing import Sequence

from .hapcut2_mec_solver import AlleleMatrix, MECSolver



def parse_cli_arguments(args=None) -> AlleleMatrix:
    parser = argparse.ArgumentParser(
        prog="HapCUT2 MEC solver",
        description="A Python wrapper around HapCUT2 to solve general minimum error correction (MEC) problems.",
        epilog="",
    )
    parser.add_argument("matrix", metavar="ALLELE_MATRIX", default="", help="The allele matrix file in json or npz format, with each row representing a fragment and each column representing a variant. Alternatively, the allele matrix can be supplied as a JSON string wrapped in quotes.")
    args = parser.parse_args(args)

    fragments: Sequence[Sequence[int]] = []
    try:
        fragments = json.loads(args.matrix)
    except Exception:
        pass
    
    if fragments:
        matrix = AlleleMatrix.from_fragments(fragments)
        return matrix

    extension: str = os.path.splitext(args.matrix)[1]
    if extension == ".json":
        with open(args.matrix, 'rt') as f:
            matrix = AlleleMatrix.from_fragments(json.load(f))
        return matrix
    elif extension == ".npz":
        matrix = AlleleMatrix.from_npz(args.matrix)
        return matrix
    else:
        raise ValueError(f"Invalid file extension: {extension!r}. Expecting .json or .npz format.")


def main(args=None, file=sys.stdout):
    matrix = parse_cli_arguments(args=args)
    result = MECSolver(matrix).solve()
    print(result.to_json(), file=file, flush=True)
