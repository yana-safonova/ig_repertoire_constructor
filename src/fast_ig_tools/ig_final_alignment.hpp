#pragma once

#include <vector>
#include <cassert>
#include <algorithm>
#include <regex>

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
seqan::String<T> consensus_hamming(const std::vector<seqan::String<T>> &reads) {
    std::vector<size_t> indices(reads.size());
    std::iota(indices.begin(), indices.end(), 0);
    return consensus_hamming(reads, indices);
}


size_t find_abundance(const std::string &s) {
    std::smatch m;
    std::regex e("(abundance:)(\\d+)$");

    if (std::regex_search(s, m, e)) {
        std::string abundance = m[2];
        return atoi(abundance.c_str());
    } else {
        return 1;
    }
}
// assert(find_abundance(">51774600-A5HPB:1:1108:11496:24798_1:N:0:TCTCGCGCATAGAGGC_abundance:666") == 666);
// assert(find_abundance(">51774600-A5HPB:1:1108:11496:24798_1:N:0:TCTCGCGCATAGAGGC_abundance:ewrwer") == 1);


template<typename T>
std::vector<size_t> find_abundances(const std::vector<T> &ids) {
    std::vector<size_t> result(ids.size());
    std::regex e("(abundance:)(\\d+)$");

    SEQAN_OMP_PRAGMA(parallel for)
        for (size_t i = 0; i < ids.size(); ++i) {
            const char *sz = seqan::toCString(ids[i]);
            std::cmatch m;

            if (std::regex_search(sz, m, e)) {
                std::string abundance = m[2];
                result[i] = atoi(abundance.c_str());
            } else {
                result[i] = 1;
            }
        }

    return result;
}

// vim: ts=4:sw=4
