#!/usr/bin/env python
# coding: utf-8

from hapcut2_mec_solver import MECSolver, solve_MEC


def test_simple_matrix():
    # Two reads with no conflicts
    allele_matrix = [[0, 0, 1, 1], [1, 1, 0, 0]]
    haplotypes = solve_MEC(allele_matrix)
    assert set(haplotypes) == {(0, 0, 1, 1), (1, 1, 0, 0)}

    # The second variant is homozygous
    allele_matrix = [[0, 1, 0], [0, 1, 0], [0, 1, 1]]
    haplotypes = solve_MEC(allele_matrix)
    assert set(haplotypes) == {(0, 1, 0), (1, 0, 1)}

    # The last read contains an error at the second variant
    allele_matrix = [[0, 1, 0], [1, 0, 1], [1, 0, 1], [1, 1, 1]]
    haplotypes = solve_MEC(allele_matrix)
    assert set(haplotypes) == {(0, 1, 0), (1, 0, 1)}


def test_ambiguous_matrix():
    # Read #3 conflicts with Read #2 at Variant #2. Without additional evidence, Variant #2 is ambiguous.
    allele_matrix = [[0, 1, 0], [1, 0, 1], [1, 1, 1]]
    haplotypes = solve_MEC(allele_matrix)
    assert set(haplotypes) == {(0, 1, 0), (1, 0, 1)}


def test_reads_with_holes():
    allele_matrix = [[0, 1, 0], [1, 0, 1], [1, -1, 1]]
    haplotypes = solve_MEC(allele_matrix)
    assert set(haplotypes) == {(0, 1, 0), (1, 0, 1)}


def test_multiple_phase_blocks():
    allele_matrix = [[0, 1, -1, -1], [1, 0, -1, -1], [-1, -1, 1, 0], [-1, -1, 0, 1]]
    haplotypes = solve_MEC(allele_matrix)
    assert set(haplotypes) == {(0, 1, 0, 1), (1, 0, 1, 0)} or set(haplotypes) == {
        (0, 1, 1, 0),
        (1, 0, 0, 1),
    }


def test_cost():
    allele_matrix = [[0, 1, 0], [1, 0, 1], [1, -1, 1], [1, 1, -1]]
    result = MECSolver(allele_matrix).solve()
    assert result.cost == 1


def test_partition():
    allele_matrix = [[0, 1, 0], [1, 0, 1], [1, -1, 1], [1, 1, -1]]
    result = MECSolver(allele_matrix).solve()
    assert result.partition in [(0, 1, 1, 0), (1, 0, 0, 1), (0, 1, 1, 1), (1, 0, 0, 0)]


def test_empty_variants():
    allele_matrix = [[0, 1, -1, -1], [1, 0, -1, -1]]
    haplotypes = solve_MEC(allele_matrix)
    assert set(haplotypes) == {(0, 1, -1, -1), (1, 0, -1, -1)}


def test_empty_fragments():
    allele_matrix = [[0, 1, 1], [1, 0, 0], [-1, -1, -1], [-1, -1, -1]]
    result = MECSolver(allele_matrix).solve()
    assert set(result.haplotypes) == {(0, 1, 1), (1, 0, 0)}
    assert result.partition in ((0, 1, -1, -1), (1, 0, -1, -1))