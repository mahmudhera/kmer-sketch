
#ifndef HASH_HPP
#define HASH_HPP

#include <cstdint>
#include <cstring>

// Lightweight 64-bit hashing utilities.
// We implement MurmurHash3_x64_128 (public domain reference) and take the low 64 bits as our 64-bit hash.
// We also provide a SplitMix64 PRNG for generating reproducible coefficient sequences (used in MinHash).

namespace hashutil {

// rotate left 64
static inline uint64_t rotl64(uint64_t x, int8_t r) {
    return (x << r) | (x >> (64 - r));
}

// FMix function from MurmurHash3 reference (public domain)
static inline uint64_t fmix64(uint64_t k) {
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33;
    return k;
}

// MurmurHash3_x64_128 (public domain), returning two 64-bit halves in out[2].
static inline void MurmurHash3_x64_128(const void * key, const int len,
                                       uint32_t seed, void * out) {
    const uint8_t * data = (const uint8_t*)key;
    const int nblocks = len / 16;

    uint64_t h1 = seed;
    uint64_t h2 = seed;

    const uint64_t c1 = 0x87c37b91114253d5ULL;
    const uint64_t c2 = 0x4cf5ad432745937fULL;

    // body
    const uint64_t * blocks = (const uint64_t *)(data);

    for (int i = 0; i < nblocks; i++) {
        uint64_t k1 = blocks[i*2+0];
        uint64_t k2 = blocks[i*2+1];

        k1 *= c1; k1 = rotl64(k1,31); k1 *= c2; h1 ^= k1;
        h1 = rotl64(h1,27); h1 += h2; h1 = h1*5+0x52dce729;

        k2 *= c2; k2 = rotl64(k2,33); k2 *= c1; h2 ^= k2;
        h2 = rotl64(h2,31); h2 += h1; h2 = h2*5+0x38495ab5;
    }

    // tail
    const uint8_t * tail = (const uint8_t*)(data + nblocks*16);
    uint64_t k1 = 0;
    uint64_t k2 = 0;

    switch(len & 15) {
    case 15: k2 ^= ((uint64_t)tail[14]) << 48;
    case 14: k2 ^= ((uint64_t)tail[13]) << 40;
    case 13: k2 ^= ((uint64_t)tail[12]) << 32;
    case 12: k2 ^= ((uint64_t)tail[11]) << 24;
    case 11: k2 ^= ((uint64_t)tail[10]) << 16;
    case 10: k2 ^= ((uint64_t)tail[ 9]) << 8;
    case  9: k2 ^= ((uint64_t)tail[ 8]) << 0;
             k2 *= c2; k2 = rotl64(k2,33); k2 *= c1; h2 ^= k2;

    case  8: k1 ^= ((uint64_t)tail[ 7]) << 56;
    case  7: k1 ^= ((uint64_t)tail[ 6]) << 48;
    case  6: k1 ^= ((uint64_t)tail[ 5]) << 40;
    case  5: k1 ^= ((uint64_t)tail[ 4]) << 32;
    case  4: k1 ^= ((uint64_t)tail[ 3]) << 24;
    case  3: k1 ^= ((uint64_t)tail[ 2]) << 16;
    case  2: k1 ^= ((uint64_t)tail[ 1]) << 8;
    case  1: k1 ^= ((uint64_t)tail[ 0]) << 0;
             k1 *= c1; k1 = rotl64(k1,31); k1 *= c2; h1 ^= k1;
    };

    // finalization
    h1 ^= len; h2 ^= len;
    h1 += h2;
    h2 += h1;

    h1 = fmix64(h1);
    h2 = fmix64(h2);

    h1 += h2;
    h2 += h1;

    ((uint64_t*)out)[0] = h1;
    ((uint64_t*)out)[1] = h2;
}

// Return 64-bit hash as the low 64 bits of MurmurHash3_x64_128.
static inline uint64_t murmurhash64(const void* data, size_t len, uint64_t seed) {
    uint64_t out[2];
    MurmurHash3_x64_128(data, (int)len, (uint32_t)seed, out);
    return out[0];
}

// SplitMix64 PRNG (for coefficients generation)
struct SplitMix64 {
    uint64_t state;
    explicit SplitMix64(uint64_t seed) : state(seed) {}
    uint64_t next() {
        uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
};

// A simple 64-bit mix (for double hashing)
static inline uint64_t mix64(uint64_t x) {
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

} // namespace hashutil

#endif // HASH_HPP
