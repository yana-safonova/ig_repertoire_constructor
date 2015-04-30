#pragma once

#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include "../include/runtime_k.hpp"

namespace string_inexact_match_index {

typedef runtime_k::RtSeq KMer;

class IndexSelector {
public:
    virtual bool NeedIncludeToIndex(const std::string & seq, size_t offset, size_t length) const = 0;
    virtual ~IndexSelector() {}
};

class AllIndexSelector : public IndexSelector {
public:
	virtual bool NeedIncludeToIndex(const std::string &, size_t, size_t) const {
		return true;
	}

    virtual ~AllIndexSelector() {
    }
};

struct SharedKMerString {
    bool operator == (const SharedKMerString & other) const {
       return seq_number_in_set == other.seq_number_in_set && pattern_start_position == other.pattern_start_position;
    }

    struct Hash {
        size_t operator () (const SharedKMerString & other) const {
            return other.seq_number_in_set * 997 + other.pattern_start_position;
        }
    };

    size_t seq_number_in_set;
    int pattern_start_position;
};

class SharedKMerStringSet {
public:
	SharedKMerStringSet(std::shared_ptr <IndexSelector> index_selector,
                        size_t k_mer_length,
                        size_t min_overlap,
                        size_t max_k_mer_occurrences);

    size_t Insert(const std::string & seq);

    void Find(const std::string & pattern, std::vector <SharedKMerString> & shared_kmer_strings) const;

private:
    size_t GetSequencesOverlap(const string & pattern, const SharedKMerString & shared_kmer_string) const;

    struct KMerMatch {
        size_t seq_number;
        size_t offset;
    };

    std::vector <size_t> sequence_lengths_;
    std::unordered_map <KMer, std::vector <KMerMatch>, KMer::hash> matches_by_kmer_;

    std::shared_ptr <IndexSelector> index_selector_;
    size_t k_mer_length_;
    size_t min_overlap_;
    size_t max_k_mer_occurrences_;
};

SharedKMerStringSet::SharedKMerStringSet(std::shared_ptr <IndexSelector> index_selector,
                                   size_t k_mer_length,
                                   size_t min_overlap,
                                   size_t max_k_mer_occurrences)
    : index_selector_(index_selector), k_mer_length_(k_mer_length), min_overlap_(min_overlap), max_k_mer_occurrences_(max_k_mer_occurrences) {
}

size_t SharedKMerStringSet::Insert(const std::string & seq) {
    size_t seq_number = sequence_lengths_.size();
    sequence_lengths_.push_back(seq.length());

    KMer kmer(0);
    for (size_t i = 0; i < k_mer_length_ && i < seq.length(); ++i) {
        kmer.pushBackThis(seq[i]);
    }
    for (size_t offset = 0; offset + k_mer_length_ <= seq.length(); kmer <<= seq[offset + k_mer_length_], ++offset) {
        if (index_selector_ -> NeedIncludeToIndex(seq, offset, k_mer_length_)) {
            matches_by_kmer_[kmer].push_back({seq_number, offset});
        // TODO OPTIMIZATION: size of vector can exceed max_k_mer_occurences
        }
    }
    return seq_number;
}

size_t SharedKMerStringSet::GetSequencesOverlap(const string & pattern, const SharedKMerString & shared_kmer_string) const {
    if (shared_kmer_string.pattern_start_position >= 0) {
        return std::min(sequence_lengths_[shared_kmer_string.seq_number_in_set] - shared_kmer_string.pattern_start_position, pattern.length());
    } else {
        return std::min(sequence_lengths_[shared_kmer_string.seq_number_in_set], pattern.length() + shared_kmer_string.pattern_start_position);
    }
}

void SharedKMerStringSet::Find(const std::string & pattern, std::vector <SharedKMerString> & shared_kmer_strings) const {
    KMer kmer(0);
    for (size_t i = 0; i < k_mer_length_ && i < pattern.length(); ++i) {
        kmer.pushBackThis(pattern[i]);
    }
    std::unordered_set <SharedKMerString, SharedKMerString::Hash> shared_kmer_strings_set;
    for (size_t offset = 0; offset + k_mer_length_ <= pattern.length(); kmer <<= pattern[offset + k_mer_length_], ++offset) {
        // TODO Skip some tricky Katya's check here
        auto matches_iterator = matches_by_kmer_.find(kmer);
        if (matches_iterator == matches_by_kmer_.end()) {
            continue;
        }
        const std::vector <KMerMatch> & matches = matches_iterator -> second;
        if (matches.size() > max_k_mer_occurrences_) {
            continue;
        }
        for (const KMerMatch & match : matches) {
            int pattern_start_position = (int)(match.offset - offset);
            SharedKMerString shared_kmer_string{match.seq_number, pattern_start_position};
            if (GetSequencesOverlap(pattern, shared_kmer_string) < min_overlap_)
                continue;
            shared_kmer_strings_set.insert(shared_kmer_string);
        }
    }

    shared_kmer_strings.resize(0);
    shared_kmer_strings.reserve(shared_kmer_strings_set.size());
    for (const auto & shared_kmer_string : shared_kmer_strings_set) {
    	shared_kmer_strings.push_back(shared_kmer_string);
    }
}

}
