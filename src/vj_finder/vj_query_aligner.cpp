#include <unordered_set>
#include <verify.hpp>

#include <core_utils.hpp>
#include "vj_query_aligner.hpp"

namespace vj_finder {
    void VJQueryAligner::CheckDbConsistencyFatal() {
        VERIFY_MSG(v_custom_db_.num_dbs() == j_custom_db_.num_dbs(), "Size of V gene DB (" << v_custom_db_.num_dbs() <<
                ") does not match with J gene DB (" << j_custom_db_.num_dbs() << ")");
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
        if(seqan::length(read) - end_of_v < algorithm_params_.filtering_params.min_j_segment_length)
            return seqan::Dna5String();
        return seqan::suffix(read, end_of_v + 1);
    }

    VJHits VJQueryAligner::Align(const core::Read &read) {
        using namespace algorithms;
        TRACE("VJ Aligner algorithm starts");
        // we can construct it once and use as a parameter of constructor
        CustomGermlineDbHelper v_kmer_index_helper(v_custom_db_);
        SubjectQueryKmerIndex<germline_utils::CustomGeneDatabase, seqan::Dna5String> v_kmer_index(
                v_custom_db_, algorithm_params_.aligner_params.word_size_v, v_kmer_index_helper);
        TRACE("Kmer index for V gene segment DB was constructed");
        PairwiseBlockAligner<germline_utils::CustomGeneDatabase, seqan::Dna5String> v_aligner(
                v_kmer_index, v_kmer_index_helper,
                CreateBlockAlignmentScoring<VJFinderConfig::AlgorithmParams::ScoringParams::VScoringParams>(
                algorithm_params_.scoring_params.v_scoring),
                CreateVBlockAlignerParams());
        TRACE("Computation of V hits");
        CustomDbBlockAlignmentHits v_aligns = v_aligner.Align(read.seq);
        TRACE(v_aligns.size() << " V hits were computed: ")
        for(auto it = v_aligns.begin(); it != v_aligns.end(); it++) {
            TRACE(v_custom_db_[it->second].name() << ", start: " << it->first.first_match_read_pos() <<
                    ", end: " << it->first.last_match_read_pos());
        }
        core::Read stranded_read = read;
        bool strand = true;
        if(algorithm_params_.aligner_params.fix_strand) {
            core::Read read_rc = read.ReverseComplement();
            CustomDbBlockAlignmentHits reverse_v_aligns = v_aligner.Align(read_rc.seq);
            if(v_aligns.BestScore() < reverse_v_aligns.BestScore()) {
                TRACE("Reverse complementary strand was selected");
                stranded_read = read_rc;
                v_aligns = reverse_v_aligns;
                strand = false;
            }
        }
        if(v_aligns.size() == 0 or !VAlignmentsAreConsistent(v_aligns))
            return VJHits(read);
        germline_utils::ChainType v_chain_type = IdentifyLocus(v_aligns);
        TRACE("V Locus was identified: " << v_chain_type);
        TRACE("Strand: " << strand);

        const germline_utils::ImmuneGeneDatabase& j_gene_db = j_custom_db_.GetConstDbByGeneType(
                germline_utils::ImmuneGeneType(v_chain_type, germline_utils::SegmentType::JoinSegment));
        TRACE("J database for locus " << v_chain_type << " consists of " << j_gene_db.size() << " gene segments");

        ImmuneGeneGermlineDbHelper j_kmer_index_helper(j_gene_db);
        SubjectQueryKmerIndex<germline_utils::ImmuneGeneDatabase, seqan::Dna5String> j_kmer_index(
                j_gene_db, algorithm_params_.aligner_params.word_size_j, j_kmer_index_helper);
        PairwiseBlockAligner<germline_utils::ImmuneGeneDatabase, seqan::Dna5String> j_aligner(
                j_kmer_index, j_kmer_index_helper,
                CreateBlockAlignmentScoring<VJFinderConfig::AlgorithmParams::ScoringParams::JScoringParams>(
                        algorithm_params_.scoring_params.j_scoring),
                CreateJBlockAlignerParams());
        auto dj_read_suffix = DefineReadJSuffix(v_aligns, stranded_read.seq);
        if(seqan::length(dj_read_suffix) == 0)
            return VJHits(read);

        TRACE("Computation of J hits");
        auto j_aligns = j_aligner.Align(dj_read_suffix);
        //for(auto it = j_aligns.begin(); it != j_aligns.end(); it++)
        //    it->first.add_read_shift(int(stranded_read.length() - seqan::length(dj_read_suffix)));
        TRACE(j_aligns.size() << " J hits were computed: ")
        for(auto it = j_aligns.begin(); it != j_aligns.end(); it++) {
            TRACE(j_gene_db[it->second].name() << ", Q start: " << it->first.first_match_read_pos() <<
            ", Q end: " << it->first.last_match_read_pos() << ", S start: " << it->first.first_match_subject_pos() <<
            ", S end: " << it->first.last_match_subject_pos());
        }
        read_archive_.UpdateReadByIndex(read.id, stranded_read.seq);
        VJHits vj_hits(read);
        for(auto it = v_aligns.begin(); it != v_aligns.end(); it++)
            vj_hits.AddVHit(VGeneHit(read, v_custom_db_[it->second], it->first, strand));
        for(auto it = j_aligns.begin(); it != j_aligns.end(); it++) {
            JGeneHit j_hit(read, j_gene_db[it->second], it->first, strand);
            j_hit.AddShift(int(read.length() - seqan::length(dj_read_suffix)));
            vj_hits.AddJHit(j_hit);
        }
        return vj_hits;
    }
}
