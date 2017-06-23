#pragma once

#include "vj_finder_config.hpp"
#include <germline_utils/germline_databases/custom_gene_database.hpp>
#include "vj_alignment_structs.hpp"
#include "block_alignment/pairwise_block_aligner.hpp"

using namespace algorithms;

namespace vj_finder {
    class VJQueryAligner {
        const VJFinderConfig::AlgorithmParams & algorithm_params_;

        core::ReadArchive &read_archive_;
        const germline_utils::CustomGeneDatabase &v_custom_db_;
        const germline_utils::CustomGeneDatabase &j_custom_db_;


        std::shared_ptr<PairwiseBlockAligner<germline_utils::CustomGeneDatabase, seqan::Dna5String> > __v_aligner = nullptr;
        
        std::vector<germline_utils::ChainType> __j_chains;
        std::vector<std::shared_ptr<PairwiseBlockAligner<germline_utils::ImmuneGeneDatabase, seqan::Dna5String> > > __j_aligners;

        void CheckDbConsistencyFatal();

        template<typename ConfigStruct>
        algorithms::BlockAlignmentScoringScheme CreateBlockAlignmentScoring(const ConfigStruct& cfg) const {
            algorithms::BlockAlignmentScoringScheme scoring;
            scoring.gap_extention_cost = cfg.gap_extention_cost;
            scoring.gap_opening_cost = cfg.gap_opening_cost;
            scoring.match_reward = cfg.match_reward;
            scoring.max_global_gap = cfg.max_global_gap;
            scoring.max_local_deletions = cfg.max_local_deletions;
            scoring.max_local_insertions = cfg.max_local_insertions;
            scoring.mismatch_extention_cost = cfg.mismatch_extention_cost;
            scoring.mismatch_opening_cost = cfg.mismatch_opening_cost;
            return scoring;
        }

        algorithms::BlockAlignerParams CreateVBlockAlignerParams() const {
                return algorithms::BlockAlignerParams(algorithm_params_.aligner_params.min_k_coverage_v,
                                                      algorithm_params_.aligner_params.max_candidates_v);
        }

        algorithms::BlockAlignerParams CreateJBlockAlignerParams() const {
            return algorithms::BlockAlignerParams(algorithm_params_.aligner_params.min_k_coverage_j,
                                                  algorithm_params_.aligner_params.max_candidates_j);
        }

        typedef algorithms::BlockAlignmentHits<germline_utils::CustomGeneDatabase> CustomDbBlockAlignmentHits;

        typedef algorithms::BlockAlignmentHits<germline_utils::ImmuneGeneDatabase> ImmuneDbBlockAlignmentHits;

        bool VAlignmentsAreConsistent(const CustomDbBlockAlignmentHits& v_alignments) const;

        germline_utils::ChainType IdentifyLocus(const CustomDbBlockAlignmentHits& v_alignments) const;

        seqan::Dna5String DefineReadJSuffix(const CustomDbBlockAlignmentHits& v_alignments,
                                            seqan::Dna5String read) const;

        template<typename SubjectDatabase, typename StringType>
        std::shared_ptr<algorithms::PairwiseBlockAligner<SubjectDatabase, StringType> > get_aligner(const algorithms::SubjectQueryKmerIndex<SubjectDatabase, StringType> &kmer_index,
                                                                                                    algorithms::KmerIndexHelper<SubjectDatabase, StringType> &kmer_index_helper,
                                                                                                    algorithms::BlockAlignmentScoringScheme scoring,
                                                                                                    algorithms::BlockAlignerParams params);

        std::shared_ptr<PairwiseBlockAligner<germline_utils::CustomGeneDatabase, seqan::Dna5String> > get_v_aligner() {
            if (!__v_aligner) {
                CustomGermlineDbHelper& v_kmer_index_helper = *(new CustomGermlineDbHelper(v_custom_db_)); //memory leak
                SubjectQueryKmerIndex<germline_utils::CustomGeneDatabase, seqan::Dna5String>& v_kmer_index =  //memory leak
                    *(new SubjectQueryKmerIndex<germline_utils::CustomGeneDatabase, seqan::Dna5String>(
                    v_custom_db_, algorithm_params_.aligner_params.word_size_v, v_kmer_index_helper));
                TRACE("Kmer index for V gene segment DB was constructed");

                __v_aligner = get_aligner(
                    v_kmer_index, v_kmer_index_helper,
                    CreateBlockAlignmentScoring<VJFinderConfig::AlgorithmParams::ScoringParams::VScoringParams>(algorithm_params_.scoring_params.v_scoring),
                    CreateVBlockAlignerParams());
            }
            return __v_aligner;
        };
        std::shared_ptr<PairwiseBlockAligner<germline_utils::ImmuneGeneDatabase, seqan::Dna5String> > get_j_aligner(germline_utils::ChainType v_chain_type) {
            auto it = find(__j_chains.begin(), __j_chains.end(), v_chain_type);
            if (it != __j_chains.end()) {
                size_t i = it - __j_chains.begin();
                return __j_aligners[i]; 
            }
            const germline_utils::ImmuneGeneDatabase& j_gene_db = j_custom_db_.GetConstDbByGeneType(
                germline_utils::ImmuneGeneType(v_chain_type, germline_utils::SegmentType::JoinSegment));

            ImmuneGeneGermlineDbHelper& j_kmer_index_helper = *(new ImmuneGeneGermlineDbHelper(j_gene_db)); //memory leak
            SubjectQueryKmerIndex<germline_utils::ImmuneGeneDatabase, seqan::Dna5String>& j_kmer_index = 
                *(new SubjectQueryKmerIndex<germline_utils::ImmuneGeneDatabase, seqan::Dna5String>( //memory leak
                    j_gene_db, algorithm_params_.aligner_params.word_size_j, j_kmer_index_helper));
            
            std::shared_ptr<PairwiseBlockAligner<germline_utils::ImmuneGeneDatabase, seqan::Dna5String> > j_aligner = get_aligner(
                j_kmer_index, j_kmer_index_helper,
                CreateBlockAlignmentScoring<VJFinderConfig::AlgorithmParams::ScoringParams::JScoringParams>(
                        algorithm_params_.scoring_params.j_scoring),
                CreateJBlockAlignerParams());

            __j_chains.push_back(v_chain_type);
            __j_aligners.push_back(j_aligner);
        
            return j_aligner;
        }
        
    public:
        VJQueryAligner(const VJFinderConfig::AlgorithmParams &algorithm_params,
                       core::ReadArchive &read_archive,
                       const germline_utils::CustomGeneDatabase &v_custom_db,
                       const germline_utils::CustomGeneDatabase &j_custom_db) :
                algorithm_params_(algorithm_params),
                read_archive_(read_archive),
                v_custom_db_(v_custom_db),
                j_custom_db_(j_custom_db) {
            CheckDbConsistencyFatal();
        }

        VJHits Align(const core::Read& read);

    private:
        DECL_LOGGER("VJQueryAligner");
    };
}
