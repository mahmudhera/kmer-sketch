import json
import os

def test_if_same():
    # read sourmash sketch hashes. this is a json file
    sourmash_sketch = json.load(open("results/fmh_sketch_sourmash", "r"))
    sourmash_hashes = sourmash_sketch[0]['signatures'][0]['mins']
    sourmash_hashes_set = set(sourmash_hashes)

    # read our implementation sketch hashes
    with open("results/fmh_sketch_our_impl", "r") as f:
        lines = f.readlines()
        our_hashes = [int(line.strip()) for line in lines if not line.startswith("#")]
    our_hashes_set = set(our_hashes)

    return sourmash_hashes_set == our_hashes_set


if __name__ == "__main__":
    # all fasta files in data/
    directory = "data/"

    # list all fasta files in directory
    fasta_files = [f for f in os.listdir(directory) if f.endswith(".fasta") or f.endswith(".fa")]

    # also list all fastq files in directory
    fasta_files += [f for f in os.listdir(directory) if f.endswith(".fastq") or f.endswith(".fq")]

    for fasta_file in fasta_files:
        fasta_file = os.path.join(directory, fasta_file)
        print(f"Testing file: {fasta_file}")

        # create fmh sketch using sourmash
        cmd = f"sourmash sketch dna {fasta_file} -o results/fmh_sketch_sourmash -f"
        print(f"Running command: {cmd}")
        os.system(cmd)

        # create fmh sketch using our implementation
        cmd = f"bin/sketch --input {fasta_file} --kmer 31 --algo fracminhash --scale 0.001 --seed 42 --output results/fmh_sketch_our_impl --canonical"
        print(f"Running command: {cmd}")
        os.system(cmd)

        # test if the two sketches are the same
        assert test_if_same(), f"Sketches do not match for {fasta_file}"
        print(f"Sketches match for {fasta_file}")

    print("All tests passed!")