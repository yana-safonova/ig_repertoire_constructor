#pragma once

#include <vector>
#include <cassert>
#include <algorithm>

#include <seqan/seq_io.h>
#include <seqan/align.h>
#include <seqan/graph_msa.h>


template<typename T = seqan::Dna5>
seqan::String<T> consensus(const std::vector<seqan::String<T>> &reads,
                           const std::vector<size_t> &indices) {
    using namespace seqan;
    using _String = seqan::String<T>;
    Align<_String> align;

    resize(rows(align), indices.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        assignSource(row(align, i), reads[indices[i]]);
    }

    globalMsaAlignment(align, SimpleScore(5, -3, -1, -3));

    String<ProfileChar<T>> profile;
    resize(profile, length(row(align, 0)));
    for (size_t rowNo = 0; rowNo < indices.size(); ++rowNo)
        for (size_t i = 0; i < length(row(align, rowNo)); ++i)
            profile[i].count[ordValue(row(align, rowNo)[i])] += 1;

    // call consensus from this string
    _String consensus;
    for (size_t i = 0; i < length(profile); ++i) {
        size_t idx = getMaxIndex(profile[i]);
        if (idx < ValueSize<T>::VALUE) {  // is not gap  TODO Check it!!
            appendValue(consensus, T(getMaxIndex(profile[i])));
        }
    }

    return consensus;
}


template<typename T = seqan::Dna5>
seqan::String<T> consensus(const std::vector<seqan::String<T>> &reads) {
    std::vector<size_t> indices(reads.size());
    std::iota(indices.begin(), indices.end(), 0);
    return consensus(reads, indices);
}


template<typename T = seqan::Dna5>
seqan::String<T> consensus_hamming(const std::vector<seqan::String<T>> &reads,
                                   const std::vector<size_t> &indices) {
    using namespace seqan;
    using _String = seqan::String<T>;

    String<ProfileChar<T>> profile;

    size_t len = 0;
    for (size_t i : indices) {
        len = std::max(len, length(reads[i]));
    }

    resize(profile, len);

    for (size_t i : indices) {
        const auto &read = reads[i];
        for (size_t j = 0; j < length(read); ++j) {
            profile[j].count[ordValue(read[j])] += 1;
        }
    }

    // call consensus from this string
    _String consensus;
    for (size_t i = 0; i < length(profile); ++i) {
        size_t idx = getMaxIndex(profile[i]);
        if (idx < ValueSize<T>::VALUE) {  // is not gap  TODO Check it!!
            appendValue(consensus, T(getMaxIndex(profile[i])));
        }
    }

    return consensus;
}


template<typename T = seqan::Dna5>
seqan::String<T> consensus_hamming_limited_coverage(const std::vector<seqan::String<T>> &reads,
                                                    const std::vector<size_t> &indices,
                                                    const std::vector<size_t> &abundances,
                                                    size_t coverage_limit = 5) {
    using namespace seqan;
    using std::vector;
    using std::max;
    using std::min;
    using _String = seqan::String<T>;

    // Compute result length
    size_t len = 0;
    for (size_t i : indices) {
        len = max(len, length(reads[i]));
    }

    String<ProfileChar<T>> profile;
    resize(profile, len);
    vector<size_t> coverage(len);

    for (size_t i : indices) {
        const auto &read = reads[i];
        for (size_t j = 0; j < length(read); ++j) {
            profile[j].count[ordValue(read[j])] += static_cast<unsigned int>(abundances[i]);
            coverage[j] += abundances[i];
        }
    }

    coverage_limit = min(coverage_limit, coverage[0]);

    _String consensus;
    for (size_t i = 0; i < len; ++i) {
        if (coverage[i] < coverage_limit) break;

        size_t idx = getMaxIndex(profile[i]);
        if (idx < ValueSize<T>::VALUE) {  // is not gap  TODO Check it!!
            appendValue(consensus, T(idx));
        }
    }

    return consensus;
}


template<typename T = seqan::Dna5>
seqan::String<T> consensus_hamming(const std::vector<seqan::String<T>> &reads) {
    std::vector<size_t> indices(reads.size());
    std::iota(indices.begin(), indices.end(), 0);
    return consensus_hamming(reads, indices);
}


template<typename T>
std::vector<size_t> find_abundances(const std::vector<T> &ids) {
    std::vector<size_t> result(ids.size());
    const std::string pat = "___size___";

    SEQAN_OMP_PRAGMA(parallel for)
    for (size_t i = 0; i < ids.size(); ++i) {
        std::string s = seqan::toCString(ids[i]);
        size_t pos = s.rfind(pat);

        if (pos == std::string::npos) {
            result[i] = 1;
        } else {
            result[i] = atoi(s.c_str() + pos + pat.length());
        }
    }

    return result;
}

// vim: ts=4:sw=4
