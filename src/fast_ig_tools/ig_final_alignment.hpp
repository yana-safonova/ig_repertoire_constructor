#define SEQAN_HAS_ZLIB 1 // TODO Set it using cmake
#define SEQAN_HAS_BZLIB 1


#include <vector>
#include <cassert>
#include <algorithm>

#include <seqan/seq_io.h>
#include <seqan/align.h>
#include <seqan/graph_msa.h>


template<typename T = seqan::Dna5>
seqan::String<T> consensus(const std::vector<seqan::String<T>> &reads,
                           const std::vector<size_t> indices) {
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
