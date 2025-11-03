#include <iostream>
#include <string>

#include "KmerScanner.hpp"
#include "FastxReader.hpp"

using namespace std;

void test_kmers(const std::string& input_path, size_t k, bool canonical = false, uint64_t seed = 42) {
    KmerScanOptions opt;
    opt.k = k;
    opt.canonical = canonical;
    opt.skip_ambiguous = true;
    opt.seed = seed;

    FastxReader reader(input_path);
    std::string header, seq;

    size_t total = 0;
    while (reader.next_record(header, seq)) {
        std::cout << ">" << header << "\n";
        scan_kmers_in_sequence(seq, opt, [&](uint64_t h) {
            std::string kmer = seq.substr(total, k);  // optional, only if you want actual substring
            std::cout << kmer << "\t" << h << "\n";
            ++total;
        });
    }

    std::cerr << "Total kmers processed: " << total << std::endl;
}

int main() {
    std::string fasta_filename = "data/my_fasta.fasta";
    int ksize = 29;
    bool canonical = false;
    test_kmers(fasta_filename, ksize, canonical);
    return 0;
}