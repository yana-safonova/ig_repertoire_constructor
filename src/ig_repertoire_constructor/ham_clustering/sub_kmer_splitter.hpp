#pragma once

#include "verify.hpp"
#include <vector>
#include <string>
#include <algorithm>
#include "io/mmapped_reader.hpp"

#include "sub_kmer_data.hpp"

#ifdef USE_GLIBCXX_PARALLEL
#include <parallel/algorithm>
#endif

namespace ham_clustering {

template <class KMerData>
class SubKMerSplitter {
    const std::string ifname_;
    const std::string ofname_;

public:
    SubKMerSplitter(const std::string &ifname, const std::string &ofname)
            : ifname_(ifname), ofname_(ofname) {}

    std::pair<size_t, size_t> split();

private:
    template<class Writer>
    void serialize(Writer &os,
                   const typename std::vector <SubKMerData <KMerData> >::iterator &start,
                   const typename std::vector <SubKMerData <KMerData> >::iterator &end) {
        size_t sz = end - start;

        os.write((char*)&sz, sizeof(sz));
        for (auto I = start, E = end; I != E; ++I)
            binary_write(os, *I);
    }

    template<class Reader>
    void deserialize(std::vector<SubKMerData<KMerData> > &res, Reader &is) {
    	res.clear();
    #if 0
    	res.shrink_to_fit();
    #else
    	std::vector<SubKMerData<KMerData> >().swap(res);
    #endif

    	size_t sz;
    	is.read((char*)&sz, sizeof(sz));
    	res.resize(sz);

    	for (size_t i = 0, e = sz; i != e; ++i)
    		binary_read(is, res[i]);
    }
};

template <class KMerData>
struct SubKMerComparator {
    bool operator()(const SubKMerData <KMerData> &lhs, const SubKMerData <KMerData> &rhs) {
        return typename KMerData::SubKMer::less2_fast()(lhs.data, rhs.data);
    }
};

template <class KMerData>
inline std::pair <size_t, size_t> SubKMerSplitter<KMerData>::split() {
    std::vector <SubKMerData <KMerData> > data;

    MMappedReader ifs(ifname_, /* unlink */ true);
    std::ofstream ofs(ofname_, std::ios::out | std::ios::binary);
    VERIFY(ofs.good());
    size_t icnt = 0, ocnt = 0;
    while (ifs.good()) {
        SubKMerComparator <KMerData> comp;

        deserialize(data, ifs);

#ifdef USE_GLIBCXX_PARALLEL
        // Explicitly force a call to parallel sort routine.
        __gnu_parallel::sort(data.begin(), data.end(), comp);
#else
        std::sort(data.begin(), data.end(), comp);
#endif
        for (auto start = data.begin(), end = data.end(); start != end;) {
            auto chunk_end = std::upper_bound(start + 1, data.end(), *start, comp);
            serialize(ofs, start, chunk_end);
            start = chunk_end;
            ocnt += 1;
        }
        icnt += 1;
    }
    VERIFY(!ofs.fail());

    ofs.close();

    return std::make_pair(icnt, ocnt);
}

}
