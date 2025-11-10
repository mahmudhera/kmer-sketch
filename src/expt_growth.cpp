#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "FastxReader.hpp"
#include "KmerScanner.hpp"
#include "Sketches.hpp"

// ----------- stable seeded 64-bit hashing of strings -----------
static inline uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
static inline uint64_t hash_string_seeded(const std::string& s, uint64_t seed) {
    uint64_t h = std::hash<std::string>{}(s);
    return splitmix64(h ^ (seed + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2)));
}

// ----------------------------- CLI helpers -----------------------------------
static void usage() {
    std::cerr
        << "Usage:\n"
        << "  mgs_similarity_experiment "
        << "--t FLOAT --metric {jaccard|cosine} "
        << "[--k K=50] [--seeds N=50] [--base_n N=1000] [--steps S=10] "
        << "[--growth {x2|x10}=x2] [--out PATH=results/mgs_similarity_experiment] "
        << "[--seed SEED=42] [--w W=64] "
        << "[--algo {maxgeom|alphamaxgeom|fracminhash|minhash|bottomk}=maxgeom] "
        << "[--alpha A=0.5] [--scale S=0.1] [--num-perm M=128]\n";
}

static std::string get_arg(std::vector<std::string>& args, const std::string& key, const std::string& def="") {
    for (size_t i=0;i+1<args.size();++i) if (args[i]==key) return args[i+1];
    return def;
}
static bool has_flag(const std::vector<std::string>& args, const std::string& key) {
    return std::find(args.begin(), args.end(), key) != args.end();
}

// ----------------------- Random string generation ----------------------------
static std::vector<std::string> generate_random_strings(size_t n, size_t length, std::mt19937_64& rng) {
    static const char alphabet[] =
        "abcdefghijklmnopqrstuvwxyz"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "0123456789";
    static const size_t A = sizeof(alphabet) - 1;

    std::uniform_int_distribution<size_t> pick(0, A - 1);
    std::vector<std::string> out;
    out.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        std::string s;
        s.resize(length);
        for (size_t j = 0; j < length; ++j) s[j] = alphabet[pick(rng)];
        out.emplace_back(std::move(s));
    }
    return out;
}

// ---------------------- Set synthesis (targets) ------------------------------
static std::pair<std::unordered_set<std::string>, std::unordered_set<std::string>>
synthesize_sets_jaccard(double t, size_t n, const std::vector<std::string>& pool, std::mt19937_64& rng) {
    size_t x = static_cast<size_t>(std::llround((2.0 * n * t) / (1.0 + t)));
    if (x > n) x = n;

    std::vector<size_t> idx(pool.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::shuffle(idx.begin(), idx.end(), rng);

    std::unordered_set<std::string> shared;
    shared.reserve(x);
    for (size_t i = 0; i < x; ++i) shared.insert(pool[idx[i]]);

    std::unordered_set<std::string> used = shared;
    std::unordered_set<std::string> A = shared, B = shared;

    size_t needA = n - x, needB = n - x;
    size_t cur = x;
    while (needA && cur < idx.size()) { const std::string& s = pool[idx[cur++]]; if (used.insert(s).second) { A.insert(s); --needA; } }
    while (needB && cur < idx.size()) { const std::string& s = pool[idx[cur++]]; if (used.insert(s).second) { B.insert(s); --needB; } }
    return {std::move(A), std::move(B)};
}

static std::pair<std::unordered_set<std::string>, std::unordered_set<std::string>>
synthesize_sets_cosine(double t, size_t n, const std::vector<std::string>& pool, std::mt19937_64& rng) {
    size_t x = static_cast<size_t>(std::llround(t * n));
    if (x > n) x = n;

    std::vector<size_t> idx(pool.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::shuffle(idx.begin(), idx.end(), rng);

    std::unordered_set<std::string> shared;
    shared.reserve(x);
    for (size_t i = 0; i < x; ++i) shared.insert(pool[idx[i]]);

    std::unordered_set<std::string> used = shared;
    std::unordered_set<std::string> A = shared, B = shared;

    size_t needA = n - x, needB = n - x;
    size_t cur = x;
    while (needA && cur < idx.size()) { const std::string& s = pool[idx[cur++]]; if (used.insert(s).second) { A.insert(s); --needA; } }
    while (needB && cur < idx.size()) { const std::string& s = pool[idx[cur++]]; if (used.insert(s).second) { B.insert(s); --needB; } }
    return {std::move(A), std::move(B)};
}

// -------------------------- MSE helper ---------------------------------------
static inline double mse(const std::vector<double>& vs, double target) {
    if (vs.empty()) return std::numeric_limits<double>::quiet_NaN();
    long double acc = 0.0L;
    for (double v : vs) {
        long double d = static_cast<long double>(v) - static_cast<long double>(target);
        acc += d * d;
    }
    return static_cast<double>(acc / vs.size());
}

// -------------------------- Sample-size helpers ------------------------------
template <typename T>
static auto has_size(const T* t) -> decltype(t->size(), bool()) { return true; }
static auto has_size(...) -> bool { return false; }

template <typename T>
static auto has_num_perm(const T* t) -> decltype(t->num_perm(), bool()) { return true; }
static auto has_num_perm(...) -> bool { return false; }

template <typename T>
static auto has_buckets(const T* t) -> decltype(t->buckets(), bool()) { return true; }
static auto has_buckets(...) -> bool { return false; }

template <typename SketchT>
static double effective_sample_size(const SketchT& s) {
    if constexpr (requires { s.size(); }) {
        return static_cast<double>(s.size());
    } else if constexpr (requires { s.num_perm(); }) {
        return static_cast<double>(s.num_perm());
    } else if constexpr (requires { s.buckets(); }) {
        // sum sizes across all buckets
        size_t total = 0;
        for (const auto& kv : s.buckets()) total += kv.second.size();
        return static_cast<double>(total);
    } else {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

// -------------------------- Estimation helpers -------------------------------
enum class Metric { Jaccard, Cosine };
struct EstResult { double estimate=std::numeric_limits<double>::quiet_NaN(); double sampleA=0.0; double sampleB=0.0; };

// MaxGeom: member jaccard/cosine
static EstResult estimate_maxgeom(const std::unordered_set<std::string>& A,
                                  const std::unordered_set<std::string>& B,
                                  uint64_t seed, Metric metric,
                                  size_t k, size_t w)
{
    MaxGeomSample sA(k, w, seed), sB(k, w, seed);
    for (const auto& x: A) sA.add_hash(hash_string_seeded(x, seed));
    for (const auto& x: B) sB.add_hash(hash_string_seeded(x, seed));
    EstResult r;
    r.estimate = (metric==Metric::Jaccard) ? sA.jaccard(sB) : sA.cosine(sB);
    r.sampleA = effective_sample_size(sA);
    r.sampleB = effective_sample_size(sB);
    return r;
}

// AlphaMaxGeom: member jaccard/cosine
static EstResult estimate_alphamaxgeom(const std::unordered_set<std::string>& A,
                                       const std::unordered_set<std::string>& B,
                                       uint64_t seed, Metric metric,
                                       double alpha, size_t w)
{
    AlphaMaxGeomSample sA(alpha, w, seed), sB(alpha, w, seed);
    for (const auto& x: A) sA.add_hash(hash_string_seeded(x, seed));
    for (const auto& x: B) sB.add_hash(hash_string_seeded(x, seed));
    EstResult r;
    r.estimate = (metric==Metric::Jaccard) ? sA.jaccard(sB) : sA.cosine(sB);
    r.sampleA = effective_sample_size(sA);
    r.sampleB = effective_sample_size(sB);
    return r;
}

// FracMinHash: static jaccard/cosine
static EstResult estimate_fracminhash(const std::unordered_set<std::string>& A,
                                      const std::unordered_set<std::string>& B,
                                      uint64_t seed, Metric metric,
                                      double scale)
{
    FracMinHash sA(scale, seed), sB(scale, seed);
    for (const auto& x: A) sA.add_hash(hash_string_seeded(x, seed));
    for (const auto& x: B) sB.add_hash(hash_string_seeded(x, seed));
    EstResult r;
    r.estimate = (metric==Metric::Jaccard) ? FracMinHash::jaccard(sA, sB)
                                           : FracMinHash::cosine(sA, sB);
    r.sampleA = effective_sample_size(sA);
    r.sampleB = effective_sample_size(sB);
    return r;
}

// MinHash: static jaccard/cosine
static EstResult estimate_minhash(const std::unordered_set<std::string>& A,
                                  const std::unordered_set<std::string>& B,
                                  uint64_t seed, Metric metric,
                                  size_t num_perm)
{
    MinHash sA(num_perm, seed), sB(num_perm, seed);
    for (const auto& x: A) sA.add_hash(hash_string_seeded(x, seed));
    for (const auto& x: B) sB.add_hash(hash_string_seeded(x, seed));
    EstResult r;
    r.estimate = (metric==Metric::Jaccard) ? MinHash::jaccard(sA, sB)
                                           : MinHash::cosine(sA, sB);
    r.sampleA = effective_sample_size(sA);
    r.sampleB = effective_sample_size(sB);
    return r;
}

// BottomK: static jaccard/cosine
static EstResult estimate_bottomk(const std::unordered_set<std::string>& A,
                                  const std::unordered_set<std::string>& B,
                                  uint64_t seed, Metric metric,
                                  size_t k)
{
    BottomK sA(k, seed), sB(k, seed);
    for (const auto& x: A) sA.add_hash(hash_string_seeded(x, seed));
    for (const auto& x: B) sB.add_hash(hash_string_seeded(x, seed));
    EstResult r;
    r.estimate = (metric==Metric::Jaccard) ? BottomK::jaccard(sA, sB)
                                           : BottomK::cosine(sA, sB);
    r.sampleA = effective_sample_size(sA);
    r.sampleB = effective_sample_size(sB);
    return r;
}

// ---------------------------- Experiment core --------------------------------
static void run_experiment(double t,
                           Metric metric,
                           const std::string& algo,
                           size_t k,
                           size_t seeds_per_size,
                           size_t base_n,
                           size_t steps,
                           const std::string& growth, // "x2" or "x10"
                           const std::string& out_path,
                           uint64_t global_seed,
                           size_t w,          // for (Alpha)MaxGeom
                           double alpha,      // AlphaMaxGeom
                           double scale,      // FracMinHash
                           size_t num_perm)   // MinHash
{
    uint64_t scale_factor = (growth=="x2") ? 2 : (growth=="x10") ? 10 : 0;
    if (!scale_factor) throw std::runtime_error("growth must be 'x2' or 'x10'");

    std::mt19937_64 rng(global_seed);

    const size_t max_size_needed = base_n * static_cast<size_t>(std::pow(scale_factor, steps ? (steps - 1) : 0)) * 2;
    const size_t pool_size = std::max<size_t>(max_size_needed, 1);
    auto universal_pool = generate_random_strings(pool_size, /*length=*/10, rng);

    std::filesystem::path outp(out_path);
    if (outp.has_parent_path() && !outp.parent_path().empty())
        std::filesystem::create_directories(outp.parent_path());

    {   // header
        std::ofstream f(out_path, std::ios::trunc);
        if (!f) throw std::runtime_error("Cannot open output file: " + out_path);
        f << "metric\tk\tstep\t|A|\t|B|\tmean_sample_size_A\tmean_sample_size_B\ttrue_sim\tmean_est\tmse\n";
    }

    std::vector<uint64_t> seeds(seeds_per_size);
    std::iota(seeds.begin(), seeds.end(), 0ULL);

    for (size_t step = 0; step < steps; ++step) {
        const size_t n = static_cast<size_t>(static_cast<long double>(base_n) * std::pow(static_cast<long double>(scale_factor), static_cast<long double>(step)));

        std::unordered_set<std::string> A, B;
        if (metric == Metric::Jaccard) {
            std::tie(A, B) = synthesize_sets_jaccard(t, n, universal_pool, rng);
        } else {
            std::tie(A, B) = synthesize_sets_cosine(t, n, universal_pool, rng);
        }

        // true similarity
        size_t inter = 0;
        if (A.size() < B.size()) {
            for (auto& s : A) if (B.count(s)) ++inter;
        }
        else {
            for (auto& s : B) if (A.count(s)) ++inter;
        }

        const double true_j = double(inter) / double(A.size() + B.size() - inter);
        const double true_c = double(inter) / std::sqrt(double(A.size()) * double(B.size()));
        const double true_sim = (metric == Metric::Jaccard) ? true_j : true_c;

        std::vector<double> ests; ests.reserve(seeds.size());
        long double sum_sA = 0.0L, sum_sB = 0.0L;

        for (uint64_t s : seeds) {
            std::cout << "Step " << step << ", n=" << n << ", seed=" << s << "\n";
            EstResult r;
            if (algo == "maxgeom") {
                r = estimate_maxgeom(A, B, s, metric, k, w);
            } else if (algo == "alphamaxgeom") {
                r = estimate_alphamaxgeom(A, B, s, metric, alpha, w);
            } else if (algo == "fracminhash") {
                r = estimate_fracminhash(A, B, s, metric, scale);
            } else if (algo == "minhash") {
                r = estimate_minhash(A, B, s, metric, num_perm);
            } else if (algo == "bottomk") {
                r = estimate_bottomk(A, B, s, metric, k);
            } else {
                throw std::runtime_error("Unsupported --algo: " + algo + " (valid: maxgeom, alphamaxgeom, fracminhash, minhash, bottomk)");
            }
            ests.push_back(r.estimate);
            sum_sA += r.sampleA;
            sum_sB += r.sampleB;
        }

        const double mean_est = ests.empty() ? std::numeric_limits<double>::quiet_NaN()
                                             : std::accumulate(ests.begin(), ests.end(), 0.0) / double(ests.size());
        const double mean_sA = double(sum_sA / seeds.size());
        const double mean_sB = double(sum_sB / seeds.size());
        const double err = mse(ests, true_sim);

        // NOTE: we keep the column name 'k' for compatibility:
        // - maxgeom/bottomk: this is K
        // - minhash: this is num_perm
        // - alphamaxgeom/fracminhash: we write 0 here (see CLI print for alpha/scale)
        size_t k_col = 0;
        if (algo == "maxgeom" || algo == "bottomk") k_col = k;
        else if (algo == "minhash") k_col = num_perm;
        else k_col = 0;

        std::ofstream f(out_path, std::ios::app);
        if (!f) throw std::runtime_error("Cannot open output file for append: " + out_path);
        f << (metric == Metric::Jaccard ? "jaccard" : "cosine") << '\t'
          << k_col << '\t'
          << step << '\t'
          << A.size() << '\t'
          << B.size() << '\t'
          << std::fixed << std::setprecision(6) << mean_sA << '\t'
          << std::fixed << std::setprecision(6) << mean_sB << '\t'
          << std::fixed << std::setprecision(6) << true_sim << '\t'
          << std::fixed << std::setprecision(6) << mean_est << '\t'
          << std::scientific << std::setprecision(6) << err
          << '\n';

        std::cerr << "Step " << step << " done: n=" << n
                  << " true=" << true_sim
                  << " mean_est=" << mean_est
                  << " mse=" << err << "\n";
    }

    std::cout << "\nResults written to:\n" << out_path << "\n";
}

// --------------------------------- main --------------------------------------
int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 1; }
    std::vector<std::string> args(argv+1, argv+argc);
    if (has_flag(args, "-h") || has_flag(args, "--help")) { usage(); return 0; }

    const std::string t_str = get_arg(args, "--t");
    const std::string metric_str = get_arg(args, "--metric");
    if (t_str.empty() || metric_str.empty()) {
        usage(); std::cerr << "\nMissing required --t or --metric.\n"; return 2;
    }

    double t = 0.0;
    try { t = std::stod(t_str); } catch (...) { std::cerr << "Invalid --t\n"; return 2; }

    Metric metric;
    if (metric_str == "jaccard" || metric_str == "Jaccard") metric = Metric::Jaccard;
    else if (metric_str == "cosine" || metric_str == "Cosine") metric = Metric::Cosine;
    else { std::cerr << "Invalid --metric (use jaccard or cosine)\n"; return 2; }

    size_t k       = static_cast<size_t>(std::stoull(get_arg(args, "--k", "50")));
    size_t seeds   = static_cast<size_t>(std::stoull(get_arg(args, "--seeds", "50")));
    size_t base_n  = static_cast<size_t>(std::stoull(get_arg(args, "--base_n", "1000")));
    size_t steps   = static_cast<size_t>(std::stoull(get_arg(args, "--steps", "10")));
    std::string growth = get_arg(args, "--growth", "x2");
    std::string out = get_arg(args, "--out", "results/mgs_similarity_experiment");
    uint64_t seed  = static_cast<uint64_t>(std::stoull(get_arg(args, "--seed", "42")));
    size_t w       = static_cast<size_t>(std::stoull(get_arg(args, "--w", "64")));
    std::string algo = get_arg(args, "--algo", "maxgeom");

    double alpha   = std::stod(get_arg(args, "--alpha", "0.5"));
    double scale   = std::stod(get_arg(args, "--scale", "0.1"));
    size_t num_perm = static_cast<size_t>(std::stoull(get_arg(args, "--num-perm", "128")));

    std::cout << "Running with the following parameters:\n";
    std::cout << "  t: " << t << "\n";
    std::cout << "  metric: " << (metric == Metric::Jaccard ? "jaccard" : "cosine") << "\n";
    std::cout << "  algo: " << algo << "\n";
    std::cout << "  k: " << k << "  (used by: maxgeom, bottomk; for minhash this is ignored in favor of --num-perm)\n";
    std::cout << "  seeds: " << seeds << "\n";
    std::cout << "  base_n: " << base_n << "\n";
    std::cout << "  steps: " << steps << "\n";
    std::cout << "  growth: " << growth << "\n";
    std::cout << "  out: " << out << "\n";
    std::cout << "  seed: " << seed << "\n";
    std::cout << "  w: " << w << "      (used by: maxgeom, alphamaxgeom)\n";
    std::cout << "  alpha: " << alpha << "  (used by: alphamaxgeom)\n";
    std::cout << "  scale: " << scale << "  (used by: fracminhash)\n";
    std::cout << "  num-perm: " << num_perm << "  (used by: minhash)\n";

    try {
        run_experiment(t, metric, algo, k, seeds, base_n, steps, growth, out, seed, w, alpha, scale, num_perm);
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 3;
    }
    return 0;
}
