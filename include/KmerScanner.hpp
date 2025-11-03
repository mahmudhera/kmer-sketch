
#ifndef KMER_SCANNER_HPP
#define KMER_SCANNER_HPP

#include <cstdint>
#include <string>
#include <functional>
#include "Util.hpp"
#include "Hash.hpp"

// KmerScanner: iterate over all k-mers of sequences and feed hashes to a callback.
// Options: skip kmers with ambiguous bases, canonicalize by reverse complement, hash seed.

struct KmerScanOptions {
    size_t k = 31;
    bool skip_ambiguous = true;
    bool canonical = false;
    uint64_t seed = 42;
};

// Calls cb(h) for each k-mer hash h encountered.
template <typename Callback>
void scan_kmers_in_sequence(const std::string& seq_raw, const KmerScanOptions& opt, Callback cb) {
    if (seq_raw.size() < opt.k) return;
    // We avoid repeated allocations by sliding a window; for canonicalization we need substrings anyway.
    const size_t n = seq_raw.size();
    std::string seq = seq_raw; // make a copy because we might uppercase in-place
    // uppercase for consistent hashing
    for (char& c : seq) {
        if (c >= 'a' && c <= 'z') c = char(c - 'a' + 'A');
    }
    for (size_t i = 0; i + opt.k <= n; ++i) {
        const char* p = &seq[i];
        bool ok = true;
        if (opt.skip_ambiguous) {
            for (size_t j = 0; j < opt.k; ++j) {
                char c = p[j];
                if (!(c=='A' || c=='C' || c=='G' || c=='T')) { ok = false; break; }
            }
        }
        if (!ok) continue;
        if (!opt.canonical) {
            uint64_t h = hashutil::murmurhash64(p, opt.k, opt.seed);
            cb(h);
        } else {
            // canonical by reverse complement: choose lexicographically smaller of kmer and its RC.
            // Build rc on the fly
            std::string rc; rc.resize(opt.k);
            for (size_t j = 0; j < opt.k; ++j) {
                char c = p[opt.k - 1 - j];
                switch (c) {
                    case 'A': rc[j] = 'T'; break;
                    case 'C': rc[j] = 'G'; break;
                    case 'G': rc[j] = 'C'; break;
                    case 'T': rc[j] = 'A'; break;
                    default: rc[j] = 'N'; ok = false; break;
                }
            }
            if (!ok && opt.skip_ambiguous) continue;
            // compare lexicographically
            bool take_rc = false;
            for (size_t j = 0; j < opt.k; ++j) {
                char a = p[j], b = rc[j];
                if (a < b) { take_rc = false; break; }
                if (a > b) { take_rc = true; break; }
            }
            if (take_rc) {
                uint64_t h = hashutil::murmurhash64(rc.data(), opt.k, opt.seed);
                cb(h);
            } else {
                uint64_t h = hashutil::murmurhash64(p, opt.k, opt.seed);
                cb(h);
            }
        }
    }
}

#endif // KMER_SCANNER_HPP
