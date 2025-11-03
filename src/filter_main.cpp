
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>

#include "SketchIO.hpp"

static void usage() {
    std::cerr << "Usage:\n"
              << "  filter --query Q.sketch --refs R1.sketch [R2.sketch ...] --metric {jaccard|cosine} --threshold X [--output OUT.tsv]\n";
}

static std::string get_arg(std::vector<std::string>& args, const std::string& key, const std::string& def="") {
    for (size_t i=0;i<args.size();++i) if (args[i]==key && i+1<args.size()) return args[i+1];
    return def;
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 1; }
    std::vector<std::string> args(argv+1, argv+argc);
    std::string qpath = get_arg(args, "--query");
    std::string metric = get_arg(args, "--metric", "jaccard");
    std::string outpath = get_arg(args, "--output", "");
    double threshold = std::stod(get_arg(args, "--threshold", "0.0"));

    // Collect all refs after '--refs' or any extra args not starting with '--'
    std::vector<std::string> refs;
    bool after_refs = false;
    for (size_t i=0;i<args.size();++i) {
        if (args[i] == "--refs") { after_refs = true; continue; }
        if (after_refs) {
            if (!args[i].empty() && args[i].rfind("--",0)!=0) refs.push_back(args[i]);
            else after_refs = false;
        }
    }
    // Also append any bare args that look like filenames
    for (const auto& s: args) if (!s.empty() && s.rfind("--",0)!=0 && s != qpath) refs.push_back(s);
    // de-dup
    std::sort(refs.begin(), refs.end()); refs.erase(std::unique(refs.begin(), refs.end()), refs.end());

    if (qpath.empty() || refs.empty()) { usage(); return 1; }

    // Load query
    VariantSketch q = load_sketch(qpath);

    // Load refs and compute scores
    struct Row { std::string path; double score; size_t inter; size_t size1; size_t size2; size_t uni; };
    std::vector<Row> rows;

    for (const auto& rp : refs) {
        if (rp == qpath) continue;
        try {
            VariantSketch r = load_sketch(rp);
            std::string why;
            if (!compatible(q, r, why)) {
                std::cerr << "Skipping incompatible reference '" << rp << "': " << why << "\n";
                continue;
            }
            double s = 0.0; size_t inter=0, uni=0, size1=0, size2=0;
            if (q.maxgeom && r.maxgeom) {
                s = (metric=="cosine") ? q.maxgeom->cosine(*r.maxgeom) : q.maxgeom->jaccard(*r.maxgeom);
                // Gather extra stats
                // For simplicity we provide sketch sizes (sum of bucket sizes).
                for (auto& kv: q.maxgeom->buckets()) size1 += kv.second.size();
                for (auto& kv: r.maxgeom->buckets()) size2 += kv.second.size();
            } else if (q.alphamaxgeom && r.alphamaxgeom) {
                s = (metric=="cosine") ? q.alphamaxgeom->cosine(*r.alphamaxgeom) : q.alphamaxgeom->jaccard(*r.alphamaxgeom);
                for (auto& kv: q.alphamaxgeom->buckets()) size1 += kv.second.size();
                for (auto& kv: r.alphamaxgeom->buckets()) size2 += kv.second.size();
            } else if (q.fracmh && r.fracmh) {
                // Jaccard and cosine on hashed sets
                if (metric=="cosine") s = FracMinHash::cosine(*q.fracmh, *r.fracmh);
                else s = FracMinHash::jaccard(*q.fracmh, *r.fracmh);
                size1 = q.fracmh->hashes().size(); size2 = r.fracmh->hashes().size();
                // compute inter/uni
                for (auto x: q.fracmh->hashes()) if (r.fracmh->hashes().count(x)) ++inter;
                uni = size1 + size2 - inter;
            } else if (q.minhash && r.minhash) {
                s = (metric=="cosine") ? MinHash::cosine(*q.minhash,*r.minhash) : MinHash::jaccard(*q.minhash,*r.minhash);
                size1 = q.minhash->num_perm(); size2 = r.minhash->num_perm();
                // counts equal components as "intersection"
                for (size_t i=0;i<q.minhash->num_perm();++i) if (q.minhash->mins()[i]==r.minhash->mins()[i]) ++inter;
                uni = q.minhash->num_perm(); // components
            } else if (q.bottomk && r.bottomk) {
                s = (metric=="cosine") ? BottomK::cosine(*q.bottomk, *r.bottomk) : BottomK::jaccard(*q.bottomk, *r.bottomk);
                size1 = q.bottomk->hashes().size(); size2 = r.bottomk->hashes().size();
                for (auto x: q.bottomk->hashes()) if (r.bottomk->hashes().count(x)) ++inter;
                uni = size1 + size2 - inter;
            } else {
                std::cerr << "Skipping '" << rp << "': internal type mismatch\n";
                continue;
            }
            if (s >= threshold) rows.push_back({rp, s, inter, size1, size2, uni});
        } catch (const std::exception& e) {
            std::cerr << "Error processing '" << rp << "': " << e.what() << "\n";
        }
    }

    // sort by score desc
    std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b){ return a.score > b.score; });

    // Output
    std::ostream* pout = &std::cout;
    std::ofstream fout;
    if (!outpath.empty()) {
        fout.open(outpath);
        if (!fout) { std::cerr << "Cannot open --output '" << outpath << "'\n"; return 2; }
        pout = &fout;
    }
    *pout << "reference\t" << metric << "_score\tintersection\tsize_query\tsize_ref\tunion\n";
    for (auto& r : rows) {
        *pout << r.path << "\t" << r.score << "\t" << r.inter << "\t" << r.size1 << "\t" << r.size2 << "\t" << r.uni << "\n";
    }
    return 0;
}
