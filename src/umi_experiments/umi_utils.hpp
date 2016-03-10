#pragma once

#include <verify.hpp>
#include <seqan/seq_io.h>
#include "../ig_tools/utils/string_tools.hpp"

class Umi {
public:
    explicit Umi(const seqan::Dna5String& umi) : umi_(umi) {}
    Umi() {}

    bool operator==(const Umi &other) const { return umi_ == other.umi_; }

    seqan::Dna5String GetString() const { return umi_; }

private:
    seqan::Dna5String umi_;
};

class Read {
public:
    explicit Read(const seqan::Dna5String& read, size_t id) : read_(read), id_(id) {}

    bool operator==(const Read& other) const { return id_ == other.id_; }

    seqan::Dna5String GetSequence() const { return read_; }
    size_t GetId() const { return id_; }

private:
    const seqan::Dna5String read_;
    const size_t id_;
};

namespace std {
    template<>
    struct hash<seqan::Dna5String> {
        size_t operator()(const seqan::Dna5String& str) const {
            size_t h = 0;
            for (auto c : str) {
                h = h * 31 + seqan::ordValue(c);
            }
            return h;
        }
    };

    template<>
    struct hash<Umi> {
        size_t operator()(const Umi& umi) const {
            size_t h = hash<seqan::Dna5String>()(umi.GetString());
            return h;
        }
    };

    template<>
    struct hash<Read> {
        size_t operator()(const Read& read) const {
            size_t h = hash<size_t>()(read.GetId());
            return h;
        }
    };
}

typedef std::shared_ptr<Umi> UmiPtr;

struct UmiPtrEquals {
    bool operator()(const UmiPtr& first, const UmiPtr& second) const {
        return *first == *second;
    }
};

struct UmiPtrHash {
    size_t operator()(const UmiPtr& umiPtr) const {
        size_t h = std::hash<Umi>()(*umiPtr);
        return h;
    }
};


void extract_barcodes_from_read_ids(const std::vector<seqan::CharString>& input_ids, std::vector<seqan::Dna5String>& umis,
                                    std::vector<seqan::DnaQString>& umi_quals);

void group_reads_by_umi(const std::vector<seqan::Dna5String>& umis, std::unordered_map<Umi, std::vector<size_t> >& umi_to_reads);
