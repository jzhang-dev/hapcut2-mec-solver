from __future__ import annotations
import argparse
import sys
import os
import json
from typing import Sequence

from .hapcut2_mec_solver import AlleleMatrix, MECSolver


def parse_cli_arguments(args: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="HapCUT2 MEC solver",
        description="A Python wrapper around HapCUT2 to solve general minimum error correction (MEC) problems.",
        epilog="",
    )
    parser.add_argument(
        "matrix",
        metavar="ALLELE_MATRIX",
        default="",
        help="The allele matrix file in json or npz format, with each row representing a fragment and each column representing a variant. Alternatively, the allele matrix can be supplied as a JSON string wrapped in quotes.",
    )
    parser.add_argument(
        "--latency-wait",
        metavar="SECONDS",
        dest="latency_wait",
        default=5,
        type=int,
        help="Number of seconds to wait after HapCUT2 exits before parsing the results.",
    )
    parser.add_argument(
        "--verbose", action="store_true", help="Print debugging information to STDERR"
    )

    parsed_args = parser.parse_args(args)
    return parsed_args


def load_allele_matrix(matrix_arg: str) -> AlleleMatrix:
    fragments: Sequence[Sequence[int]] = []
    try:
        return AlleleMatrix.from_json_string(matrix_arg)
    except Exception:
        pass

    extension: str = os.path.splitext(matrix_arg)[1]
    if extension == ".json":
        return AlleleMatrix.from_json(matrix_arg)
    elif extension == ".npz":
        return AlleleMatrix.from_npz(matrix_arg)
    else:
        raise ValueError(
            f"Invalid file extension: {extension!r}. Expecting .json or .npz format."
        )


def main(args: Sequence[str] | None = None) -> None:
    parsed_args = parse_cli_arguments(args=args)
    if parsed_args.verbose:
        print("Loading allele matrix", file=sys.stderr, flush=True)
    matrix = load_allele_matrix(parsed_args.matrix)
    if parsed_args.verbose:
        n_row, n_col = matrix.shape
        print(
            f"Loaded allele matrix with {n_row} rows and {n_col} columns",
            file=sys.stderr,
            flush=True,
        )
        print(
            f"Initializing MEC solver",
            file=sys.stderr,
            flush=True,
        )
    solver = MECSolver(matrix)
    result = solver.solve(
        verbose=parsed_args.verbose, latency_wait=parsed_args.latency_wait
    )
    print(result.to_json(), file=sys.stdout, flush=True)
