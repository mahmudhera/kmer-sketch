
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
              << "[--k K] [--b B] [--w W] [--alpha A] [--scale S] [--num-perm M] [--seed SEED] "
              << "[--canonical] [--keep-ambiguous] --output OUT\n";
}



static void detailed_usage() {
    std::cerr << "Options:\n"
              << "  --input FILE            Input FASTA/FASTQ file\n"
              << "  --kmer N                K-mer size (default: 31)\n"
              << "  --algo ALGO             Sketching algorithm to use: maxgeom, alphamaxgeom, fracminhash, minhash, bottomk\n"
              << "  --k K                   (bottomk) Sketch size K (default: 1000)\n"
              << "  --b B                   (maxgeom) Bucket capacity B (default: 90)\n"
              << "  --w W                   (maxgeom, alphamaxgeom) Max number of buckets W (default: 64)\n"
              << "  --alpha A               (alphamaxgeom) Alpha parameter (default: 0.45)\n"
              << "  --scale S               (fracminhash) Scale parameter (default: 0.001)\n"
              << "  --num-perm K            (minhash) Number of permutations K (default: 1000)\n"
              << "  --seed SEED             Random seed (default: 42)\n"
              << "  --canonical             Use canonical k-mers\n"
              << "  --keep-ambiguous        Keep ambiguous k-mers (default: skip them)\n"
              << "  --output OUT            Output sketch file\n"
              << "  --help, -h              Show this help message\n";
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

    // check if --help is requested
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h") {
            usage();
            return 0;
        } else if (arg == "--help") {
            detailed_usage();
            return 0;
        }
    }

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
        size_t B = (size_t)std::stoull(get_arg(args, "--b", "90"));
        size_t W = (size_t)std::stoull(get_arg(args, "--w", "64"));
        MaxGeomSample sketch(B, W, seed);
        // scan
        FastxReader reader(inpath);
        std::string header, seq;
        while (reader.next_record(header, seq)) {
            scan_kmers_in_sequence(seq, opt, [&](uint64_t h){ sketch.add_hash(h); });
        }
        sketch.write(out, kmer);
    } else if (algo == "alphamaxgeom") {
        double alpha = std::stod(get_arg(args, "--alpha", "0.45"));
        size_t W = (size_t)std::stoull(get_arg(args, "--w", "64"));
        AlphaMaxGeomSample sketch(alpha, W, seed);
        FastxReader reader(inpath);
        std::string header, seq;
        while (reader.next_record(header, seq)) {
            scan_kmers_in_sequence(seq, opt, [&](uint64_t h){ sketch.add_hash(h); });
        }
        sketch.write(out, kmer);
    } else if (algo == "fracminhash") {
        double scale = std::stod(get_arg(args, "--scale", "0.001"));
        FracMinHash sketch(scale, seed);
        FastxReader reader(inpath);
        std::string header, seq;
        while (reader.next_record(header, seq)) {
            scan_kmers_in_sequence(seq, opt, [&](uint64_t h){ sketch.add_hash(h); });
        }
        sketch.write(out, kmer);
    } else if (algo == "minhash") {
        size_t M = (size_t)std::stoull(get_arg(args, "--num-perm", "1000"));
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
