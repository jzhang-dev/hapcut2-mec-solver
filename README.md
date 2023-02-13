# HapCUT2 MEC solver

Solve general minimum error correction (MEC) problems using [HapCUT2](https://github.com/vibansal/HapCUT2).

## Usage

HapCUT2 MEC solver uses the [JSON format](https://en.wikipedia.org/wiki/JSON) for input and output. The input JSON file should contain a two-dimentional matrix, with each row representing a fragment, and each column representing a variant. Allowed values are `0` (Allele 0), `1` (Allele 1) and `-1` (unknown allele). 

For example, the following is a valid JSON input file with three variants and four fragments:

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
# The input JSON can be supplied via STDIN
docker run -i --rm jzhang0246/hapcut2-mec-solver:latest python -m hapcut2_mec_solver < example.json
# Alternatively, the input JSON can be used directly as an argment surrounded by single or double quotes. 
docker run -i --rm jzhang0246/hapcut2-mec-solver:latest python -m hapcut2_mec_solver '[[1, 0, 1], [0, 1, 0]]'
```

Run tests:

```sh
docker run -i --rm jzhang0246/hapcut2-mec-solver:latest pytest
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

>>> allele_matrix = [[0, 0, 1, 1], [1, 1, 0, 0]]
>>> haplotypes = solve_MEC(allele_matrix)
>>> assert set(haplotypes) == {(0, 0, 1, 1), (1, 1, 0, 0)}
True

>>> allele_matrix = [[0, 1, 0], [1, 0, 1], [1, -1, 1], [1, 1, -1]]
>>> result = MECSolver(allele_matrix).solve()
>>> assert result.cost == 1
True
>>> assert result.partition in [(0, 1, 1, 0), (1, 0, 0, 1), (0, 1, 1, 1), (1, 0, 0, 0)]
True
```

Run tests:

```sh
cd hapcut2-mec-solver
pytest
```


## Notes

- HapCUT2 assumes all variants are heterozygous by default even if only one allele is seen in the fragments. For example, fragments `[[1, 0, 1], [1, 0, 1], [0, 1, 1]]` will be solved as haplotypes `[[1, 0, 1], [0, 1, 0]]`. Although a [`--call_homozygous` option](https://github.com/vibansal/HapCUT2/blob/66ee827f9130fa64ff93044227702804308f1650/hapcut2-src/optionparser.c#L201) is availble in HapCUT2, the prior probablity of a variant being homozygous is [hard-coded to be extremely small](https://github.com/vibansal/HapCUT2/blob/66ee827f9130fa64ff93044227702804308f1650/hapcut2-src/optionparser.c#L23), preventing a variant to be called homozygous in most cases. See [this issue](https://github.com/vibansal/HapCUT2/issues/110) for further discussions.  
- Fragments that contain only `-1` will be assigned `-1` in the output partitions. Variants that contain only `-1` will show as `-1` in the output haplotypes. 
- HapCUT2 harnesses a heuristic algorithm to appoach the $\mathcal{NP}$-hard MEC problem. HapCUT2 is not guaranteed to find the theoretically optimal solution. In contrast, [WhatsHap](https://whatshap.readthedocs.io/en/latest/) always finds the theoretically optimal solution with the lowest cost, but has $\mathit{O}(exp(n))$ time complexity, where $n$ is the maximum coverage of the variants. See also [WhatsHap wMEC solver](https://github.com/jzhang-dev/whatshap-wmec-solver).

## References

Edge, P., Bafna, V., & Bansal, V. (2017). HapCUT2: Robust and accurate haplotype assembly for diverse sequencing technologies. Genome Research, 27(5), 801–812. https://doi.org/10.1101/gr.213462.116



