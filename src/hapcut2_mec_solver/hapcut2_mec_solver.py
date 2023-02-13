#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
from typing import Sequence, IO, Iterator
import collections
import os
from dataclasses import dataclass
import json
import subprocess
import tempfile
import numpy as np
import pandas as pd
from scipy.sparse import csr_array, save_npz, load_npz


# Type aliases
Fragment = Sequence[int]
Haplotype = Sequence[int]
# AlleleMatrix = Sequence[Fragment]


class AlleleMatrix:
    def __init__(self, sparse_matrix: csr_array):
        self._sparse_matrix = sparse_matrix

    @classmethod
    def from_fragments(cls, fragments: Sequence[Fragment]):
        dense_matrix = np.array(fragments, dtype=np.int8)
        sparse_matrix: csr_array = csr_array(dense_matrix + 1, dtype=np.int8)
        return cls(sparse_matrix)

    def __getitem__(self, key) -> Sequence[int]:
        row = self._sparse_matrix.getrow(key).toarray() - 1
        return row.tolist()[0]

    def to_npz(self, file: str | IO) -> None:
        save_npz(file, self._sparse_matrix)

    @classmethod
    def from_npz(cls, file: str | IO):
        return cls(load_npz(file))

    def __len__(self) -> int:
        return self._sparse_matrix.shape[0]

    def __iter__(self) -> Iterator[Sequence[int]]:
        for i in range(len(self)):
            yield self[i]


@dataclass(eq=True)
class _MECSolverResult:
    haplotypes: tuple[Haplotype, Haplotype]
    partition: Sequence[int]
    cost: float

    def to_json(self):
        return json.dumps(
            dict(
                haplotypes=[list(haplotype) for haplotype in self.haplotypes],
                partition=list(self.partition),
                cost=self.cost,
            )
        )


class MECSolver:
    def __init__(self, matrix: AlleleMatrix):
        self.matrix: AlleleMatrix = matrix
        empty_variants: list[int] = []
        for j in range(self.n_variant):
            if all(matrix[i][j] < 0 for i in range(self.n_fragment)):
                empty_variants.append(j)
        self._empty_variants = empty_variants


    @classmethod
    def from_fragments(cls, fragments: Sequence[Fragment]):
        matrix = AlleleMatrix.from_fragments(fragments)
        return cls(matrix)


    @property
    def n_variant(self) -> int:
        return len(self.matrix[0])

    @property
    def n_fragment(self) -> int:
        return len(self.matrix)

    @staticmethod
    def _make_vcf(
        n_variant: int,
        output_path: str,
        *,
        ref_alleles: Sequence[str] | None = None,
        alt_alleles: Sequence[str] | None = None,
        genotypes: Sequence[str] | None = None,
    ) -> None:
        VCF_HEADER = [
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            "Sample_0",
        ]
        chromosome = ["ref"] * n_variant
        pos = [i * 100 + 1 for i in range(n_variant)]
        variant_id = ["."] * n_variant
        ref = ["A"] * n_variant if ref_alleles is None else ref_alleles
        alt = ["T"] * n_variant if alt_alleles is None else alt_alleles
        qual = [60] * n_variant
        filter_ = ["."] * n_variant
        info = ["."] * n_variant
        format_ = ["GT"] * n_variant
        gt = ["."] * n_variant if genotypes is None else genotypes
        data = {
            k: v
            for k, v in zip(
                VCF_HEADER,
                [
                    chromosome,
                    pos,
                    variant_id,
                    ref,
                    alt,
                    qual,
                    filter_,
                    info,
                    format_,
                    gt,
                ],
            )
        }
        df = pd.DataFrame(data)
        df.to_csv(output_path, index=False, sep="\t")

    def _make_fragments(self, output_path: str, *, fragment_names=None) -> None:
        n_fragment, n_variant = self.n_fragment, self.n_variant
        matrix = self.matrix
        _fragment_names = (
            [f"Read_{i+1}" for i in range(n_fragment)]
            if fragment_names is None
            else fragment_names
        )

        FRAGMENT_TYPE = 0  # 0 for normal, 1 for HiC
        BARCODE = -1
        MATE_INDEX = -1

        fragment_lines: list[str] = []
        for i in range(n_fragment):
            fragment: Sequence[int] = matrix[i]
            variant_indices: list[int] = [
                i for i, allele in enumerate(fragment) if allele >= 0
            ]
            block_start_indices: list[int] = [
                j for j in variant_indices if j == 0 or j - 1 not in variant_indices
            ]

            fragment_string = "".join([str(a) if a >= 0 else " " for a in fragment])
            blocks: list[str] = [
                block for block in fragment_string.strip(" ").split(" ") if block
            ]
            quality_string = "." * len(variant_indices)
            block_string = " ".join(
                [f"{s+1} {b}" for s, b in zip(block_start_indices, blocks)]
            )

            line = f"{len(blocks)} {_fragment_names[i]} {FRAGMENT_TYPE} {BARCODE} {MATE_INDEX} {block_string} {quality_string}\n"
            fragment_lines.append(line)

        with open(output_path, "w") as f:
            f.write("".join(fragment_lines))

    def _run_hapcut2(
        self,
        fragments_path: str,
        vcf_path: str,
        output_path: str,
        *,
        prune=False,
        call_homozygous=False,
    ) -> None:
        directory = os.path.commonprefix([fragments_path, vcf_path, output_path])
        command = [
            "hapcut2",
            "--fragments",
            fragments_path,
            "--VCF",
            vcf_path,
            "--output",
            output_path,
            "--outvcf",
            "0",
            "--call_homozygous",
            "1" if call_homozygous else "0",
            "--new_format",
            "1",
            "--verbose",
            "0",
            "--error_analysis_mode",
            str(int(not prune)),
        ]
        process = subprocess.run(command, capture_output=True, encoding="utf-8")
        # if process.returncode != 0:
        #     raise RuntimeError(f"Failed to run HapCUT2: \n{process.stderr}")

    def _parse_hapcut2_result(self, file_path: str) -> tuple[Haplotype, Haplotype]:
        haplotype_0: list[int] = []
        haplotype_1: list[int] = []

        with open(file_path, "r") as f:
            for line in f:
                if line.startswith("BLOCK"):
                    continue
                elif line.startswith("*"):
                    continue
                else:
                    columns = line.strip("\n").split("\t")
                    haplotype_0.append(int(columns[1]) if columns[1] != "-" else -1)
                    haplotype_1.append(int(columns[2]) if columns[2] != "-" else -1)

        for j in self._empty_variants:
            haplotype_0.insert(j, -1)
            haplotype_1.insert(j, -1)

        return (tuple(haplotype_0), tuple(haplotype_1))

    @staticmethod
    def _get_cost(haplotype: Haplotype, fragment: Fragment) -> int:
        return sum(
            fragment_allele >= 0 and haplotype_allele != fragment_allele
            for haplotype_allele, fragment_allele in zip(haplotype, fragment)
        )

    def _partition_fragments(self, haplotypes) -> tuple[Sequence[int], float]:
        n_fragment, n_variant = self.n_fragment, self.n_variant
        n_haplotype = len(haplotypes)
        partition: list[int] = []
        total_cost: float = 0
        for fragment in self.matrix:
            if all(allele < 0 for allele in fragment):
                partition.append(-1)
                continue
            min_cost = float("inf")
            haplotype_index: int = -1
            for i, haplotype in enumerate(haplotypes):
                cost = self._get_cost(haplotype, fragment)
                if cost < min_cost:
                    min_cost = cost
                    haplotype_index = i
            partition.append(haplotype_index)
            total_cost += min_cost
        return tuple(partition), total_cost

    def solve(self, *, call_homozygous=False) -> _MECSolverResult:
        with tempfile.TemporaryDirectory() as temp_directory:
            vcf_path = os.path.join(temp_directory, "variants.vcf")
            fragments_path = os.path.join(temp_directory, "fragments.txt")
            output_path = os.path.join(temp_directory, "hapcut2.txt")
            self._make_vcf(self.n_variant, vcf_path)
            self._make_fragments(fragments_path)
            self._run_hapcut2(
                fragments_path, vcf_path, output_path, call_homozygous=call_homozygous
            )
            haplotypes = self._parse_hapcut2_result(output_path)
            partition, cost = self._partition_fragments(haplotypes)
            return _MECSolverResult(
                haplotypes=haplotypes, partition=partition, cost=cost
            )


def solve_MEC(
    fragments: Sequence[Fragment], *, call_homozygous=False
) -> tuple[Sequence[int], Sequence[int]]:
    solver = MECSolver.from_fragments(fragments)
    result = solver.solve(call_homozygous=call_homozygous)
    haplotype_1, haplotype_2 = result.haplotypes
    return tuple(haplotype_1), tuple(haplotype_2)
