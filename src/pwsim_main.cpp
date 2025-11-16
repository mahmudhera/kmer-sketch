
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include "SketchIO.hpp"

static void usage() {
    std::cerr << "Usage:\n"
              << "  pwsimilarity --metric {jaccard|cosine} --output OUT.tsv SKETCH1 SKETCH2 [SKETCH3 ...]\n";
}

static void detailed_usage() {
    std::cerr << "Options:\n"
              << "  --metric METRIC         Similarity metric to compute: jaccard (default), cosine\n"
              << "  --output OUT.tsv        Output TSV file for pairwise similarities (default: pairs.tsv)\n"
              << "  SKETCH1 SKETCH2 ...     Input sketch files to compare\n";
}

static std::string get_arg(std::vector<std::string>& args, const std::string& key, const std::string& def="") {
    for (size_t i=0;i<args.size();++i) if (args[i]==key && i+1<args.size()) return args[i+1];
    return def;
}

int main(int argc, char** argv) {

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

    if (argc < 3) { usage(); return 1; }

    std::vector<std::string> args(argv+1, argv+argc);
    std::string metric = get_arg(args, "--metric", "jaccard");
    std::string outpath = get_arg(args, "--output", "pairs.tsv");

    // collect all filenames (args not starting with '--' and not values of known flags)
    std::vector<std::string> files;
    for (size_t i=0;i<args.size();++i) {
        if (args[i].rfind("--",0)==0) { ++i; continue; } // skip flag and its value (rough)
        files.push_back(args[i]);
    }
    if (files.size() < 2) { usage(); return 1; }

    // load all sketches
    std::vector<VariantSketch> sketches;
    sketches.reserve(files.size());
    for (const auto& p : files) sketches.push_back(load_sketch(p));

    // check compatibility
    for (size_t i=1;i<sketches.size();++i) {
        std::string why;
        if (!compatible(sketches[0], sketches[i], why)) {
            std::cerr << "Incompatible sketches: '" << files[0] << "' vs '" << files[i] << "': " << why << "\n";
            return 2;
        }
    }

    std::ofstream out(outpath);
    if (!out) { std::cerr << "Cannot open output file: " << outpath << "\n"; return 3; }
    out << "sketch1\tsketch2\t" << metric << "_score\n";

    // compute pairwise
    for (size_t i=0;i<sketches.size();++i) {
        for (size_t j=i+1;j<sketches.size();++j) {
            double s = 0.0; size_t inter=0, uni=0, size1=0, size2=0;
            const auto& a = sketches[i]; const auto& b = sketches[j];
            if (a.maxgeom && b.maxgeom) {
                s = (metric=="cosine") ? a.maxgeom->cosine(*b.maxgeom) : a.maxgeom->jaccard(*b.maxgeom);
                for (auto& kv: a.maxgeom->buckets()) size1 += kv.second.size();
                for (auto& kv: b.maxgeom->buckets()) size2 += kv.second.size();
            } else if (a.alphamaxgeom && b.alphamaxgeom) {
                s = (metric=="cosine") ? a.alphamaxgeom->cosine(*b.alphamaxgeom) : a.alphamaxgeom->jaccard(*b.alphamaxgeom);
                for (auto& kv: a.alphamaxgeom->buckets()) size1 += kv.second.size();
                for (auto& kv: b.alphamaxgeom->buckets()) size2 += kv.second.size();
            } else if (a.fracmh && b.fracmh) {
                if (metric=="cosine") s = FracMinHash::cosine(*a.fracmh, *b.fracmh);
                else s = FracMinHash::jaccard(*a.fracmh, *b.fracmh);
                size1 = a.fracmh->hashes().size(); size2 = b.fracmh->hashes().size();
                for (auto x: a.fracmh->hashes()) if (b.fracmh->hashes().count(x)) ++inter;
                uni = size1 + size2 - inter;
            } else if (a.minhash && b.minhash) {
                s = (metric=="cosine") ? MinHash::cosine(*a.minhash, *b.minhash) : MinHash::jaccard(*a.minhash, *b.minhash);
                size1 = a.minhash->num_perm(); size2 = b.minhash->num_perm();
                for (size_t t=0;t<a.minhash->num_perm();++t) if (a.minhash->mins()[t]==b.minhash->mins()[t]) ++inter;
                uni = a.minhash->num_perm();
            } else if (a.bottomk && b.bottomk) {
                s = (metric=="cosine") ? BottomK::cosine(*a.bottomk, *b.bottomk) : BottomK::jaccard(*a.bottomk, *b.bottomk);
                size1 = a.bottomk->hashes().size(); size2 = b.bottomk->hashes().size();
                for (auto x: a.bottomk->hashes()) if (b.bottomk->hashes().count(x)) ++inter;
                uni = size1 + size2 - inter;
            } else {
                std::cerr << "Internal type mismatch between " << files[i] << " and " << files[j] << "\n";
                continue;
            }
            out << files[i] << "\t" << files[j] << "\t" << s << "\n";
        }
    }
    std::cerr << "Wrote: " << outpath << "\n";
    return 0;
}
