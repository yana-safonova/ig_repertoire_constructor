#define SEQAN_HAS_ZLIB 1 // TODO Set it using cmake
#define SEQAN_HAS_BZLIB 1


#include <vector>
#include <cassert>
#include <algorithm>
#include <unordered_map>

#include <seqan/seq_io.h>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace seqan;

using std::cout;

template<typename T = seqan::Dna5>
seqan::String<T> consensus(const std::vector<seqan::String<T>> &reads) {
  using _String = seqan::String<T>;
  Align<_String> align;

  resize(rows(align), reads.size());
  for (size_t i = 0; i < reads.size(); ++i) {
    assignSource(row(align, i), reads[i]);
  }

  globalMsaAlignment(align, SimpleScore(5, -3, -1, -3));

  String<ProfileChar<T>> profile;
  resize(profile, length(row(align, 0)));
  for (size_t rowNo = 0; rowNo < reads.size(); ++rowNo)
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


int main()
{
  // some variangs of sonic hedgehog exon 1
  std::vector<Dna5String> reads {
    // gi|2440284|dbj|AB007129.1| Oryzias latipes
    "GCGGGTCACTGAGGGCTGGGATGAGGACGGCCACCACTTCGAGGAGTCCCTTCACTACGAGGGCAGGGCC"
      "GTGGACATCACCACGTCAGACAGGGACAAGAGCAAGTACGGCACCCTGTCCAGACTGGCGGTGGAAGCTG"
      "GGTTCGACTGGGTCTACTATGAGTCCAAAGCGCACATCCACTGCTCTGTGAAAGCAGAAAGCTCAGTCGC"
      "TGCAAAGTCGGGCGGTTGCTTCCCAGGATCCTCCACGGTCACCCTGGAAAATGGCACCCAGAGGCCCGTC"
      "AAAGATCTCCAACCCGGGGACAGAGTACTGGCCGCGGATTACGACGGAAACCCGGTTTATACCGACTTCA"
      "TCATGTTCAA",
    // gi|1731488|gb|U51350.1|DDU51350 Devario devario
    "CTACGGCAGAAGAAGACATCCGAAAAAGCTGACACCTCTCGCCTACAAGCAGTTCATACCTAATGTCGCG"
      "GAGAAGACCTTAGGGGCCAGCGGCAGATACGAGGGCAAGATAACGCGCAATTCGGAGAGATTTAAAGAAC"
      "TTACTCCAAATTACAATCCCGACATTATCTTTAAGGATGAGGAGAACACG",
    // gi|1731504|gb|U51352.1|PTU51352 Puntius tetrazona
    "CTACGGCAGAAGAAGACATCCCAAGAAGCTGACACCTCTCGCCTACAAGCAGTTTATACCTAATGTCGCG"
      "GAGAAGACCTTAGGGGCCAGCGGCAGATACGAGGGCAAGATCACGCGCAATTCGGAGAGATTTAAAGAAC"
      "TTACTCCAAATTACAATCCCGACATTATCTTTAAGGATGAGGAGAACACT",
    // gi|54399708|gb|AY642858.1| Bos taurus
    "TGCTGCTGCTGGCGAGATGTCTGCTGGTGCTGCTTGTCTCCTCGCTGTTGATGTGCTCGGGGCTGGCGTG"
      "CGGACCCGGCAGGGGATTTGGCAAGAGGCGGAACCCCAAAAAGCTGACCCCTTTAGCCTACAAGCAGTTT"
      "ATCCCCAACGTGGCGGAGAAGACCCTAGGGGCCAGTGGAAGATATGAGGGGAAGATCACCAGAAACTCAG"
      "AGCGATTTAAGGAACTCACCCCCAATTACAACCC"
  };

  auto cons = consensus(reads);

  cout << cons << std::endl;
  return 0;
}
