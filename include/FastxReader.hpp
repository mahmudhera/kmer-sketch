
#ifndef FASTX_READER_HPP
#define FASTX_READER_HPP

#include <string>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <iostream>

// Minimal FASTA/FASTQ reader with streaming of sequences.
// - FASTA: '>' header, 1+ sequence lines (may be split across multiple lines).
// - FASTQ: '@' header, next line sequence, '+' line, next line qualities.
// Non-ACGT bases are allowed; k-mer scanner may decide to skip.

class FastxReader {
public:
    explicit FastxReader(const std::string& path) : in_(path) {
        if (!in_) throw std::runtime_error("Failed to open input file: " + path);
        // Peek first char to guess format lazily.
    }

    bool next_record(std::string& header, std::string& sequence) {
        header.clear();
        sequence.clear();
        std::string line;
        // Skip until non-empty line
        char c;
        while (true) {
            if (!std::getline(in_, line)) return false;
            if (!line.empty()) break;
        }
        c = line[0];
        if (c == '>') { // FASTA
            header = line.substr(1);
            // Read until next header or EOF
            std::streampos pos = in_.tellg();
            while (std::getline(in_, line)) {
                if (!line.empty() && (line[0] == '>' || line[0] == '@')) {
                    // push back by resetting position to start of this line
                    in_.seekg(pos);
                    break;
                }
                sequence += line;
                pos = in_.tellg();
            }
            return true;
        } else if (c == '@') { // FASTQ
            header = line.substr(1);
            std::string seq, plus, qual;
            if (!std::getline(in_, seq)) throw std::runtime_error("Broken FASTQ: missing sequence after header @" + header);
            if (!std::getline(in_, plus)) throw std::runtime_error("Broken FASTQ: missing '+' line for @" + header);
            if (!std::getline(in_, qual)) throw std::runtime_error("Broken FASTQ: missing qualities for @" + header);
            if (plus.empty() || plus[0] != '+') throw std::runtime_error("Broken FASTQ: '+' line malformed for @" + header);
            sequence = seq;
            return true;
        } else {
            // Not a FASTX header; throw error
            throw std::runtime_error("Invalid FASTX format: unexpected line start '" + std::string(1, c) + "'");
        }
    }

private:
    std::ifstream in_;
};

#endif // FASTX_READER_HPP
