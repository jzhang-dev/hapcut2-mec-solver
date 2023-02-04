import argparse
import sys
import json

from .hapcut2_mec_solver import AlleleMatrix, MECSolver

def load_stdin():
    return sys.stdin.read()

def parse_cli_arguments() -> AlleleMatrix:
    parser = argparse.ArgumentParser(
        prog="HapCUT2 MEC solver",
        description="A Python wrapper around HapCUT2 to solve general minimum error correction (MEC) problems.",
        epilog="",
    )
    parser.add_argument("matrix_json", metavar="ALLELE_MATRIX", nargs='?', default="")
    args = parser.parse_args()
    matrix_json = args.matrix_json or load_stdin()
    if matrix_json is "":
        raise ValueError("Please specify ALLELE_MATRIX as an argument or through STDIN.")
    matrix: AlleleMatrix = json.loads(matrix_json)
    return matrix


def main():
    matrix = parse_cli_arguments()
    result = MECSolver(matrix).solve()
    print(result.to_json())