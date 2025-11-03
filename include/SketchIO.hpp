
#ifndef SKETCH_IO_HPP
#define SKETCH_IO_HPP

#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <memory>
#include "Sketches.hpp"

// A small loader that inspects the header to detect algorithm, then parses accordingly.
// We represent a "variant" type via simple struct with pointers. Callers must manage which pointer is non-null.

struct VariantSketch {
    std::string algo;
    size_t kmer_size = 0;
    uint64_t seed = 42;
    // One of the following will be non-null.
    std::unique_ptr<MaxGeomSample> maxgeom;
    std::unique_ptr<AlphaMaxGeomSample> alphamaxgeom;
    std::unique_ptr<FracMinHash> fracmh;
    std::unique_ptr<MinHash> minhash;
    std::unique_ptr<BottomK> bottomk;
};

inline VariantSketch load_sketch(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open sketch file: " + path);
    std::string first_header;
    std::string line;
    std::string algo;
    // Peek headers
    std::streampos start = in.tellg();
    while (std::getline(in, line)) {
        if (line.size()>2 && line[0]=='#') {
            auto pos = line.find('=');
            if (pos != std::string::npos) {
                std::string key = line.substr(2, pos-2);
                std::string val = line.substr(pos+1);
                if (key == "algo") { algo = val; break; }
            }
        } else if (!line.empty()) {
            break;
        }
    }
    // reset to beginning
    in.clear(); in.seekg(0);
    VariantSketch var;
    if (algo == "MaxGeom") {
        size_t kmer=0; uint64_t seed=0;
        auto s = MaxGeomSample::read(in, kmer, seed);
        var.algo = "MaxGeom"; var.kmer_size = kmer; var.seed = seed;
        var.maxgeom = std::make_unique<MaxGeomSample>(s);
    } else if (algo == "AlphaMaxGeom") {
        size_t kmer=0; uint64_t seed=0; double alpha=0.5;
        auto s = AlphaMaxGeomSample::read(in, kmer, seed, alpha);
        var.algo = "AlphaMaxGeom"; var.kmer_size = kmer; var.seed = seed;
        var.alphamaxgeom = std::make_unique<AlphaMaxGeomSample>(s);
    } else if (algo == "FracMinHash") {
        size_t kmer=0; uint64_t seed=0; double scale=1.0;
        auto s = FracMinHash::read(in, kmer, seed, scale);
        var.algo = "FracMinHash"; var.kmer_size = kmer; var.seed = seed;
        var.fracmh = std::make_unique<FracMinHash>(s);
    } else if (algo == "MinHash") {
        size_t kmer=0; uint64_t seed=0; size_t num_perm=0;
        auto s = MinHash::read(in, kmer, seed, num_perm);
        var.algo = "MinHash"; var.kmer_size = kmer; var.seed = seed;
        var.minhash = std::make_unique<MinHash>(s);
    } else if (algo == "BottomK") {
        size_t kmer=0; uint64_t seed=0; size_t kk=0;
        auto s = BottomK::read(in, kmer, seed, kk);
        var.algo = "BottomK"; var.kmer_size = kmer; var.seed = seed;
        var.bottomk = std::make_unique<BottomK>(s);
    } else {
        throw std::runtime_error("Unknown or missing 'algo' in sketch file: " + path);
    }
    return var;
}

// Compatibility check (algo, k-mer size, seed, and algorithm-specific params).
// Because our serialized files fix all needed params in the header, and we reconstruct the sketch from the file,
// simply checking 'algo', 'k-mer size', 'seed' is sufficient. (Algorithm-specific parameter mismatches will
// have led to different object types with different internal sizes / parameters.)
inline bool compatible(const VariantSketch& a, const VariantSketch& b, std::string& why) {
    if (a.algo != b.algo) { why = "different algo: " + a.algo + " vs " + b.algo; return false; }
    if (a.kmer_size != b.kmer_size) { why = "kmer_size mismatch"; return false; }
    if (a.seed != b.seed) { why = "hash_seed mismatch"; return false; }
    // For algorithms with additional params, we could compare sizes by re-serializing headers,
    // but since our derived classes are concrete, we check important fields explicitly:
    if (a.maxgeom && b.maxgeom) {
        if (a.maxgeom->k() != b.maxgeom->k() || a.maxgeom->w() != b.maxgeom->w()) { why = "MaxGeom params differ"; return false; }
    }
    if (a.alphamaxgeom && b.alphamaxgeom) {
        if (a.alphamaxgeom->w() != b.alphamaxgeom->w() || a.alphamaxgeom->alpha() != b.alphamaxgeom->alpha()) { why = "AlphaMaxGeom params differ"; return false; }
    }
    if (a.fracmh && b.fracmh) {
        if (a.fracmh->scale() != b.fracmh->scale()) { why = "FracMinHash scale differs"; return false; }
    }
    if (a.minhash && b.minhash) {
        if (a.minhash->num_perm() != b.minhash->num_perm()) { why = "MinHash num_perm differs"; return false; }
    }
    if (a.bottomk && b.bottomk) {
        if (a.bottomk->k() != b.bottomk->k()) { why = "BottomK k differs"; return false; }
    }
    return true;
}

#endif // SKETCH_IO_HPP
