import os
import random
import hashlib
import statistics
import math

# test/test_maxgeom.py

DATA_DIR = "data"
OUTPUT_DIR = "outputs"
FASTA_PATH = os.path.join(DATA_DIR, "synthetic_test_fasta.fasta")
KMER_SIZE = 31
NUM_SEEDS = 50


def gen_unique_kmers(n, k, seed=42):
    """Generate n unique DNA kmers of length k deterministically."""
    rnd = random.Random(seed)
    bases = ["A", "C", "G", "T"]
    seen = set()
    kmers = []
    while len(kmers) < n:
        s = "".join(rnd.choice(bases) for _ in range(k))
        if s not in seen:
            seen.add(s)
            kmers.append(s)
    return kmers


def write_fasta(path, kmers):
    with open(path, "w") as fh:
        for i, seq in enumerate(kmers):
            fh.write(f">kmer_{i}\n")
            fh.write(seq + "\n")


def run_alpha_max_geom_and_get_num_hashes(fasta_filename, seed, aMGS_alpha):
    cmd = f"bin/sketch --input {fasta_filename} --kmer 31 --algo alphamaxgeom --alpha {aMGS_alpha} --canonical --output outputs/test_alpha_max_geom_sample --seed {seed}"
    os.system(cmd)
    
    # read num of lines in output sketch file
    with open("outputs/test_alpha_max_geom_sample", "r") as f:
        lines = f.readlines()
        num_hashes = len(lines) - 8  # subtract header lines
    
    return num_hashes


def test_alpha_maxgeom_sampling(num_kmers, kmer_size, aMGS_alpha, num_seeds):
    # Generate fasta with NUM_KMERS kmers length K
    kmers = gen_unique_kmers(num_kmers, kmer_size, seed=42)
    write_fasta(FASTA_PATH, kmers)
    assert os.path.exists(FASTA_PATH)

    # Run maxgeom multiple times and collect sample sizes
    sample_sizes = []
    for seed in range(num_seeds):
        num_hashes = run_alpha_max_geom_and_get_num_hashes(FASTA_PATH, seed, aMGS_alpha)
        sample_sizes.append(num_hashes)

    assert sample_sizes, "No samples were collected"

    # Compute and print statistics
    mean_size = statistics.mean(sample_sizes)
    stdev_size = statistics.stdev(sample_sizes)
    return mean_size, stdev_size



if __name__ == "__main__":
    num_kmers_list = [10000, 20000, 30000, 40000, 50000]
    aMGS_alpha_list = [0.25, 0.3, 0.4, 0.5]
    for num_kmers in num_kmers_list:
        for aMGS_alpha in aMGS_alpha_list:
            r = test_alpha_maxgeom_sampling(num_kmers, KMER_SIZE, aMGS_alpha, NUM_SEEDS)
            print(f"num_kmers: {num_kmers}, aMGS_alpha: {aMGS_alpha} => mean sample size: {r[0]:.2f}, stdev: {r[1]:.2f}")