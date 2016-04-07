#include <unordered_set>
#include <verify.hpp>

#include <core_utils.hpp>
#include "vj_query_aligner.hpp"

namespace vj_finder {
    void VJQueryAligner::CheckDbConsistencyFatal() {
        VERIFY_MSG(v_custom_db_.size() == j_custom_db_.size(), "Size of V gene DB (" << v_custom_db_.size() <<
                                                               ") does not match with J gene DB (" <<
                                                               j_custom_db_.size() << ")");
    }

    bool VJQueryAligner::VAlignmentsAreConsistent(const VJQueryAligner::CustomDbBlockAlignmentHits &v_alignments) const {
        std::unordered_set<germline_utils::ImmuneGeneType, germline_utils::ImmuneGeneTypeHasher> gene_types;
        for(size_t i = 0; i < v_alignments.size(); i++)
            gene_types.insert(v_custom_db_[v_alignments[i].second].GeneType());
        return gene_types.size() == 1;
    }

    germline_utils::ChainType VJQueryAligner::IdentifyLocus(const CustomDbBlockAlignmentHits &v_alignments) const {
        if(v_alignments.size() == 0)
            return germline_utils::ChainType();
        size_t index0 = v_alignments[0].second;
        return v_custom_db_[index0].Chain();
    }

    seqan::Dna5String VJQueryAligner::DefineReadJSuffix(const CustomDbBlockAlignmentHits &v_alignments,
                                                        seqan::Dna5String read) const {
        size_t end_of_v = core::max_map(v_alignments.cbegin(), v_alignments.cend(),
                               [](const CustomDbBlockAlignmentHits::IndicedPairwiseBlockAlignment &align) ->
                                       size_t { return align.first.last_match_read_pos(); });
        return seqan::suffix(read, end_of_v + 1);
    }

    VJHits VJQueryAligner::Align(const core::Read &read) {
        using namespace algorithms;
        // we can construct it once and use as a parameter of constructor
        CustomGermlineDbHelper v_kmer_index_helper(v_custom_db_);
        SubjectQueryKmerIndex<germline_utils::CustomGeneDatabase, seqan::Dna5String> v_kmer_index(
                v_custom_db_, algorithm_params_.aligner_params.word_size_v, v_kmer_index_helper);
        PairwiseBlockAligner<germline_utils::CustomGeneDatabase, seqan::Dna5String> v_aligner(
                v_kmer_index, v_kmer_index_helper,
                CreateBlockAlignmentScoring<vjf_config::AlgorithmParams::ScoringParams::VScoringParams>(
                algorithm_params_.scoring_params.v_scoring),
                CreateVBlockAlignerParams());
        CustomDbBlockAlignmentHits v_aligns = v_aligner.Align(read.seq);
        core::Read stranded_read = read;
        bool strand = true;
        if(algorithm_params_.query_params.fix_strand) {
            core::Read read_rc = read; // todo: implement RC
            CustomDbBlockAlignmentHits reverse_v_aligns = v_aligner.Align(read_rc.seq);
            if(v_aligns.BestScore() < reverse_v_aligns.BestScore()) {
                stranded_read = read_rc;
                v_aligns = reverse_v_aligns;
            }
        }
        if(v_aligns.size() == 0 or !VAlignmentsAreConsistent(v_aligns))
            return VJHits(read);
        germline_utils::ChainType v_chain_type = IdentifyLocus(v_aligns);
        const germline_utils::ImmuneGeneDatabase& j_gene_db = j_custom_db_.GetDbByGeneType(
                germline_utils::ImmuneGeneType(v_chain_type, germline_utils::SegmentType::JoinSegment));

        ImmuneGeneGermlineDbHelper j_kmer_index_helper(j_gene_db);
        SubjectQueryKmerIndex<germline_utils::ImmuneGeneDatabase, seqan::Dna5String> j_kmer_index(
                j_gene_db, algorithm_params_.aligner_params.word_size_j, j_kmer_index_helper);
        PairwiseBlockAligner<germline_utils::ImmuneGeneDatabase, seqan::Dna5String> j_aligner(
                j_kmer_index, j_kmer_index_helper,
                CreateBlockAlignmentScoring<vjf_config::AlgorithmParams::ScoringParams::JScoringParams>(
                        algorithm_params_.scoring_params.j_scoring),
                CreateJBlockAlignerParams());
        auto j_aligns = j_aligner.Align(DefineReadJSuffix(v_aligns, stranded_read.seq));
        VJHits vj_hits(read);
        for(auto it = v_aligns.begin(); it != v_aligns.end(); it++)
            vj_hits.AddVHit(VGeneHit(read, v_custom_db_[it->second], it->first, strand));
        for(auto it = j_aligns.begin(); it != j_aligns.end(); it++)
            vj_hits.AddJHit(JGeneHit(read, j_gene_db[it->second], it->first, strand));
        return vj_hits;
    }
}
