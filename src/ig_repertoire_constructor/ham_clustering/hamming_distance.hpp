#pragma once

namespace ham_clustering {

template <class KMerData>
unsigned GetHammingDistance(const typename KMerData::KMer & kmer1,
                            const typename KMerData::KMer & kmer2,
                            unsigned max_allowed_dist) {
    VERIFY(kmer1.size() == kmer2.size());

    unsigned dist = 0;
    for (size_t i = 0; i < kmer1.size(); ++i) {
        if (kmer1[i] != kmer2[i]) {
            ++dist;
            if (dist > max_allowed_dist) {
                break;
            }
        }
    }
    return dist;
}

}
