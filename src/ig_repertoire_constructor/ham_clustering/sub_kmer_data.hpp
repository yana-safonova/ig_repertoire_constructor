#pragma once

namespace ham_clustering {

template <class KMerData>
struct SubKMerData {
	uint64_t idx;
    typename KMerData::SubKMer data;
};

//static_assert(sizeof(SubKMerData) == 16, "Too big SubKMer");    TODO ???

template <class Reader, class KMerData>
inline void binary_read(Reader &is, SubKMerData <KMerData> &s) {
	is.read((char*)&s.idx, sizeof(s.idx));
	KMerData::binary_read_subkmer(is, s.data);
}

template <class Writer, class KMerData>
inline void binary_write(Writer &os, const SubKMerData <KMerData> &s) {
    os.write((char*)&s.idx, sizeof(s.idx));
    KMerData::binary_write_subkmer(os, s.data);
}

}
