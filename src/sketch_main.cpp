
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <unordered_map>

#include "FastxReader.hpp"
#include "KmerScanner.hpp"
#include "Sketches.hpp"

static void usage() {
    std::cerr << "Usage:\n"
              << "  sketch --input FILE --kmer N --algo {maxgeom|alphamaxgeom|fracminhash|minhash|bottomk} "
              << "[--k K] [--w W] [--alpha A] [--scale S] [--num-perm M] [--seed SEED] "
              << "[--canonical] [--keep-ambiguous] --output OUT\n";
}

// Simple arg parsing
static std::string get_arg(std::vector<std::string>& args, const std::string& key, const std::string& def="") {
    for (size_t i=0;i<args.size();++i) if (args[i]==key && i+1<args.size()) return args[i+1];
    return def;
}
static bool has_flag(const std::vector<std::string>& args, const std::string& key) {
    for (size_t i=0;i<args.size();++i) if (args[i]==key) return true;
    return false;
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 1; }
    std::vector<std::string> args(argv+1, argv+argc);
    std::string inpath = get_arg(args, "--input");
    std::string algo = get_arg(args, "--algo");
    std::string outpath = get_arg(args, "--output");
    if (inpath.empty() || algo.empty() || outpath.empty()) { usage(); return 1; }
    size_t kmer = (size_t)std::stoull(get_arg(args, "--kmer", "31"));
    uint64_t seed = (uint64_t)std::stoull(get_arg(args, "--seed", "42"));
    bool canonical = has_flag(args, "--canonical");
    bool keep_ambiguous = has_flag(args, "--keep-ambiguous"); // default: skip ambiguity

    KmerScanOptions opt;
    opt.k = kmer; opt.seed = seed; opt.canonical = canonical; opt.skip_ambiguous = !keep_ambiguous;

    // Instantiate sketch
    std::ofstream out(outpath);
    if (!out) { std::cerr << "Cannot open output file: " << outpath << "\n"; return 2; }

    if (algo == "maxgeom") {
        size_t K = (size_t)std::stoull(get_arg(args, "--k", "64"));
        size_t W = (size_t)std::stoull(get_arg(args, "--w", "64"));
        MaxGeomSample sketch(K, W, seed);
        // scan
        FastxReader reader(inpath);
        std::string header, seq;
        while (reader.next_record(header, seq)) {
            scan_kmers_in_sequence(seq, opt, [&](uint64_t h){ sketch.add_hash(h); });
        }
        sketch.write(out, kmer);
    } else if (algo == "alphamaxgeom") {
        double alpha = std::stod(get_arg(args, "--alpha", "0.5"));
        size_t W = (size_t)std::stoull(get_arg(args, "--w", "64"));
        AlphaMaxGeomSample sketch(alpha, W, seed);
        FastxReader reader(inpath);
        std::string header, seq;
        while (reader.next_record(header, seq)) {
            scan_kmers_in_sequence(seq, opt, [&](uint64_t h){ sketch.add_hash(h); });
        }
        sketch.write(out, kmer);
    } else if (algo == "fracminhash") {
        double scale = std::stod(get_arg(args, "--scale", "0.1"));
        FracMinHash sketch(scale, seed);
        FastxReader reader(inpath);
        std::string header, seq;
        while (reader.next_record(header, seq)) {
            scan_kmers_in_sequence(seq, opt, [&](uint64_t h){ sketch.add_hash(h); });
        }
        sketch.write(out, kmer);
    } else if (algo == "minhash") {
        size_t M = (size_t)std::stoull(get_arg(args, "--num-perm", "128"));
        MinHash sketch(M, seed);
        FastxReader reader(inpath);
        std::string header, seq;
        while (reader.next_record(header, seq)) {
            scan_kmers_in_sequence(seq, opt, [&](uint64_t h){ sketch.add_hash(h); });
        }
        sketch.write(out, kmer);
    } else if (algo == "bottomk") {
        size_t K = (size_t)std::stoull(get_arg(args, "--k", "1000"));
        BottomK sketch(K, seed);
        FastxReader reader(inpath);
        std::string header, seq;
        while (reader.next_record(header, seq)) {
            scan_kmers_in_sequence(seq, opt, [&](uint64_t h){ sketch.add_hash(h); });
        }
        sketch.write(out, kmer);
    } else {
        std::cerr << "Unknown --algo: " << algo << "\n";
        usage(); return 1;
    }

    return 0;
}
