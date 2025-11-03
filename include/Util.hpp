
#ifndef UTIL_HPP
#define UTIL_HPP

#include <string>
#include <algorithm>
#include <cctype>
#include <vector>

namespace util {

inline std::string to_upper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::toupper(c); });
    return s;
}

inline bool is_acgt(char c) {
    switch (c) {
        case 'A': case 'C': case 'G': case 'T':
        case 'a': case 'c': case 'g': case 't':
            return true;
        default: return false;
    }
}

inline char comp_base(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'a': return 't';
        case 'c': return 'g';
        case 'g': return 'c';
        case 't': return 'a';
        default: return 'N';
    }
}

inline std::string revcomp(const std::string& s) {
    std::string rc; rc.reserve(s.size());
    for (auto it = s.rbegin(); it != s.rend(); ++it) rc.push_back(comp_base(*it));
    return rc;
}

inline bool contains_only_acgt(const std::string& s) {
    for (char c : s) if (!is_acgt(c)) return false;
    return true;
}

// Trim whitespace (both ends)
inline std::string trim(const std::string& s) {
    size_t start = 0;
    while (start < s.size() && std::isspace((unsigned char)s[start])) ++start;
    size_t end = s.size();
    while (end > start && std::isspace((unsigned char)s[end-1])) --end;
    return s.substr(start, end - start);
}

} // namespace util

#endif // UTIL_HPP
