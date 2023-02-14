# HapCUT2 MEC solver

Solve general minimum error correction (MEC) problems using [HapCUT2](https://github.com/vibansal/HapCUT2).

## Usage

### Input and output

The input for HapCUT2 MEC solver is a *allele matrix*, which is a two-dimensional matrix, with each row representing a *fragment*, and each column representing a *locus*. 

For an allele matrix with *n* rows and *m* columns, the output includes:
- a *haplotype* matrix with 2 rows and *m* columns;
- a *partitioning* vector of length *n*;
- and the minimum *cost*.



The input allele matrix can be supplied using either the [JSON](https://en.wikipedia.org/wiki/JSON) format, or as a [.npz file](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.load_npz.html#scipy-sparse-load-npz) containing a [sparse CSR matrix](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html#scipy-sparse-csr-matrix). In both cases, the output will be a specified using the JSON format. 

For JSON input, the allowed values are `0` (Allele 0), `1` (Allele 1) and `-1` (unknown allele), following the convention of MEC problems. For sparse matrix input, the allowed values are `1` (Allele 0), `2` (Allele 1) and `0` (unknown allele), following the convention of sparse matrices. The output always uses `{-1, 0, 1}` values for alleles.

For example, the following is a valid JSON input file with three loci and four fragments:

```json
[
    [ 0,  1,  0], 
    [ 1,  0,  1], 
    [ 1, -1,  1],
    [ 1,  1, -1]
]
```

The output for the input matrix above will be:

```json
{"haplotypes": [[0, 1, 0], [1, 0, 1]], "partition": [0, 1, 1, 0], "cost": 1}
```


### Docker

Pull image from Docker Hub:

```sh
docker pull jzhang0246/hapcut2-mec-solver:latest
```

Example usage: 

```sh
# The input allele matrix can be supplied as a JSON file
docker run \
  --rm \
  -v "$(pwd)"/tests:/data \
  jzhang0246/hapcut2-mec-solver:latest \
    python -m hapcut2_mec_solver \
      /data/test.json

# ... or as a npz file containing a sparse matrix 
docker run \
  --rm \
  -v "$(pwd)"/tests:/data \
  jzhang0246/hapcut2-mec-solver:latest \
    python -m hapcut2_mec_solver \
      /data/test.npz


# Alternatively, supply the matrix as a JSON string surrounded by single or double quotes. 
docker run \
  --rm 
  jzhang0246/hapcut2-mec-solver:latest \
    python -m hapcut2_mec_solver \
      '[[1, 0, 1], [0, 1, 0]]'
```

Run tests:

```sh
docker run --rm jzhang0246/hapcut2-mec-solver:latest pytest
```

### Python integration

If HapCUT2 has already been installed and is available under `$PATH` as `hapcut2`, then HapCUT2 MEC solver can be used as a Python library.

Install using `pip`:

```sh
pip install git+https://github.com/jzhang-dev/hapcut2-mec-solver
```

Example usage:

```py
>>> from hapcut2_mec_solver import MECSolver, solve_MEC

>>> fragments = [[0, 0, 1, 1], [1, 1, 0, 0]]
>>> haplotypes = solve_MEC(fragments)
>>> set(haplotypes) == {(0, 0, 1, 1), (1, 1, 0, 0)}
True

>>> fragments = [[0, 1, 0], [1, 0, 1], [1, -1, 1], [1, 1, -1]]
>>> result = MECSolver.from_fragments(fragments).solve()
>>> result.cost == 1
True
>>> result.partition in [(0, 1, 1, 0), (1, 0, 0, 1), (0, 1, 1, 1), (1, 0, 0, 0)]
True
```

Run tests:

```sh
git clone https://github.com/jzhang-dev/hapcut2-mec-solver
cd hapcut2-mec-solver
pip install -v .
pip install pytest
pytest
```


## Notes

- HapCUT2 assumes all loci are heterozygous by default even if only one allele is seen in the fragments. For example, fragments `[[1, 0, 1], [1, 0, 1], [0, 1, 1]]` will be solved as haplotypes `[[1, 0, 1], [0, 1, 0]]`. Although a [`--call_homozygous`](https://github.com/vibansal/HapCUT2/blob/66ee827f9130fa64ff93044227702804308f1650/hapcut2-src/optionparser.c#L201) option is availble in HapCUT2, the prior probablity of a locus being homozygous is [hard-coded to be extremely small](https://github.com/vibansal/HapCUT2/blob/66ee827f9130fa64ff93044227702804308f1650/hapcut2-src/optionparser.c#L23), preventing a locus to be called homozygous in most cases. See [this issue](https://github.com/vibansal/HapCUT2/issues/110) for further discussions.  
- Fragments that contain only `-1` will be assigned `-1` in the output partitions. Loci that contain only `-1` will show as `-1` in the output haplotypes. 
- HapCUT2 harnesses a heuristic algorithm to appoach the $\mathcal{NP}$-hard MEC problem. HapCUT2 is not guaranteed to find the theoretically optimal solution. In contrast, [WhatsHap](https://whatshap.readthedocs.io/en/latest/) always finds the theoretically optimal solution with the lowest cost, but has $\mathit{O}(\exp(n))$ time complexity, where $n$ is the maximum coverage of the loci. See also [WhatsHap wMEC solver](https://github.com/jzhang-dev/whatshap-wmec-solver).
- See [tests](https://github.com/jzhang-dev/hapcut2-mec-solver/tree/main/tests) for more examples about the behavior of the HapCUT2 algorithm. 

## References

Edge, P., Bafna, V., & Bansal, V. (2017). HapCUT2: Robust and accurate haplotype assembly for diverse sequencing technologies. Genome Research, 27(5), 801–812. https://doi.org/10.1101/gr.213462.116



