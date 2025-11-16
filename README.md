# MaxGeomHash & α-MaxGeomHash  
A C++ implementation of MaxGeomHash and α-MaxGeomHash for 
- k-mer sketching of FASTA/FASTQ files, and 
- pairwise similarity estimation from k-mer sketches.

## Overview
This repository provides efficient C++ implementations of **MaxGeomHash** and **α-MaxGeomHash**, two hashing–based sketching algorithms. The tools support several other sketching methods -- including FracMinHash, MinHash (bottom-k). From FASTA or FASTQ files, the tools can compute k-mer sketches using any of these sketching methods. The tools also allow for rapid pairwise similarity estimation using a number of sketches.

All dependencies are either header-only or included with the repository. Compilation requires only `make`.

---

## Build Instructions
Clone the repository and build:

```bash
make
```

The executables are generated in the bin/ directory. Kindly make sure to add the bin/ directory in your PATH variable.

## Programs

Two main programs are included:
- `sketch`
- `pwsimilarity`

### `sketch`
This program creates the sketches. Following are the arguments.

| Argument            | Argument Type | Default Value | What It Means |
|---------------------|---------------|----------------|----------------|
| `--input FILE`      | string (path) | **required**   | Input sequence file in FASTA or FASTQ format. |
| `--kmer N`          | integer       | 31             | k-mer size to break sequences into. |
| `--algo ALGO`       | string        | **required**   | Sketching algorithm to use: `maxgeom`, `alphamaxgeom`, `fracminhash`, `minhash`, `bottomk`. |
| `--k K`             | integer       | 1000           | (bottom-k) Sketch size *K* for bottom-k hashing. |
| `--b B`             | integer       | 90             | (maxgeom) Bucket capacity *B* used in MaxGeomHash. |
| `--w W`             | integer       | 64             | (maxgeom, alphamaxgeom) Maximum number of buckets *W*. |
| `--alpha A`         | float         | 0.45           | (alphamaxgeom) Alpha parameter for α-MaxGeomHash. |
| `--scale S`         | float         | 0.001          | (fracminhash) Scale parameter controlling sampling probability. |
| `--num-perm K`      | integer       | 1000           | (minhash) Number of permutations for classical MinHash. |
| `--seed SEED`       | integer       | 42             | Random seed for reproducibility. |
| `--canonical`       | flag          | false          | Treat each k-mer as canonical (min of forward/reverse complement). |
| `--keep-ambiguous`  | flag          | false          | Keep k-mers containing ambiguous bases instead of skipping. |
| `--output OUT`      | string (path) | **required**   | Output sketch file path. |
| `--help`, `-h`      | flag          | —              | Display help message and exit. |

## `pwsimilarity`

This program computes pairwise similarity from a list of sketch files. The arguments are as follows:

| Argument               | Argument Type      | Default Value | What It Means |
|------------------------|--------------------|----------------|----------------|
| `--metric METRIC`      | string             | jaccard        | Similarity metric to compute between sketches. Options: `jaccard`, `cosine`. |
| `--output OUT.tsv`     | string (path)      | pairs.tsv      | Output TSV file containing pairwise similarity results. |
| `SKETCH1 SKETCH2 ...`  | list of file paths | **required**   | Input sketch files to compare pairwise. |
