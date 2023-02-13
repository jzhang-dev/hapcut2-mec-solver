import argparse
import sys
import os
import json
from typing import Sequence

from .hapcut2_mec_solver import AlleleMatrix, MECSolver


def load_stdin():
    return sys.stdin.read()


def parse_cli_arguments(args=None) -> AlleleMatrix:
    parser = argparse.ArgumentParser(
        prog="HapCUT2 MEC solver",
        description="A Python wrapper around HapCUT2 to solve general minimum error correction (MEC) problems.",
        epilog="",
    )
    parser.add_argument("matrix", metavar="ALLELE_MATRIX", nargs="?", default="")
    args = parser.parse_args(args)

    fragments: Sequence[Sequence[int]] = []
    if args.matrix == "":
        stdin = load_stdin()
        if not stdin:
            raise RuntimeError(
            "Please specify ALLELE_MATRIX as an argument or through STDIN."
        )
        fragments = json.loads(load_stdin())
        matrix = AlleleMatrix.from_fragments(fragments)
        return matrix
    
    
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


    raise RuntimeError("Invalid arguments")


def main(args=None, file=sys.stdout):
    matrix = parse_cli_arguments(args=args)
    result = MECSolver(matrix).solve()
    print(result.to_json(), file=file, flush=True)
