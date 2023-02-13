#!/usr/bin/env python
# coding: utf-8

from hapcut2_mec_solver import MECSolver, solve_MEC, AlleleMatrix

def test_readme():
    from hapcut2_mec_solver import MECSolver, solve_MEC

    fragments = [[0, 0, 1, 1], [1, 1, 0, 0]]
    haplotypes = solve_MEC(fragments)
    assert set(haplotypes) == {(0, 0, 1, 1), (1, 1, 0, 0)}

    fragments = [[0, 1, 0], [1, 0, 1], [1, -1, 1], [1, 1, -1]]
    result = MECSolver.from_fragments(fragments).solve()
    assert result.cost == 1
    assert result.partition in [(0, 1, 1, 0), (1, 0, 0, 1), (0, 1, 1, 1), (1, 0, 0, 0)]


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


def test_homozygous_variants():
    # See https://github.com/vibansal/HapCUT2/issues/110
    # See https://github.com/vibansal/HapCUT2/blob/66ee827f9130fa64ff93044227702804308f1650/hapcut2-src/optionparser.c#L23

    # The second variant is homozygous
    allele_matrix = [[0, 1, 0] for _ in range(500)] + [[1, 1, 1] for _ in range(500)]
    haplotypes = solve_MEC(allele_matrix, call_homozygous=True)
    assert set(haplotypes) == {(0, 1, 0), (1, 1, 1)}


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
    result = MECSolver.from_fragments(allele_matrix).solve()
    assert result.cost == 1


def test_partition():
    allele_matrix = [[0, 1, 0], [1, 0, 1], [1, -1, 1], [1, 1, -1]]
    result = MECSolver.from_fragments(allele_matrix).solve()
    assert result.partition in [(0, 1, 1, 0), (1, 0, 0, 1), (0, 1, 1, 1), (1, 0, 0, 0)]


def test_empty_variants():
    allele_matrix = [[0, 1, -1, -1], [1, 0, -1, -1]]
    haplotypes = solve_MEC(allele_matrix)
    assert set(haplotypes) == {(0, 1, -1, -1), (1, 0, -1, -1)}


def test_empty_fragments():
    allele_matrix = [[0, 1, 1], [1, 0, 0], [-1, -1, -1], [-1, -1, -1]]
    result = MECSolver.from_fragments(allele_matrix).solve()
    assert set(result.haplotypes) == {(0, 1, 1), (1, 0, 0)}
    assert result.partition in ((0, 1, -1, -1), (1, 0, -1, -1))



def test_sparse_matrix():
    import tempfile, os

    fragments = [[0, 0, 1, 1], [1, 1, 0, 0], [-1, -1, 0, 0]]
    matrix1 = AlleleMatrix.from_fragments(fragments)
    assert len(matrix1) == len(fragments)

    with tempfile.TemporaryDirectory() as temp_directory:
        npz_path = os.path.join(temp_directory, "test.npz")
        matrix1.to_npz(npz_path)
        matrix2 = AlleleMatrix.from_npz(npz_path)
    
    assert len(matrix2) == len(fragments)

    result1 = MECSolver(matrix1).solve()
    result2 = MECSolver(matrix2).solve()

    assert result1 == result2



