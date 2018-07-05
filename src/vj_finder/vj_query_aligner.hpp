#pragma once

#include "vj_finder_config.hpp"
#include <germline_utils/germline_databases/custom_gene_database.hpp>
#include "vj_alignment_structs.hpp"
#include "block_alignment/pairwise_block_aligner.hpp"

namespace vj_finder {
    class VJQueryAligner {

        const VJFinderConfig::AlgorithmParams & algorithm_params_;

        //core::ReadArchive &read_archive_;
        const germline_utils::CustomGeneDatabase &v_custom_db_;
        const germline_utils::CustomGeneDatabase &j_custom_db_;

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
                                                                                                    algorithms::BlockAlignerParams params) {
            return std::make_shared<algorithms::PairwiseBlockAligner<SubjectDatabase, StringType> >(kmer_index, kmer_index_helper, scoring, params);
        }

    public:
        VJQueryAligner(const VJFinderConfig::AlgorithmParams &algorithm_params,
                       //core::ReadArchive &read_archive,
                       const germline_utils::CustomGeneDatabase &v_custom_db,
                       const germline_utils::CustomGeneDatabase &j_custom_db) :
                algorithm_params_(algorithm_params),
                //read_archive_(read_archive),
                v_custom_db_(v_custom_db),
                j_custom_db_(j_custom_db) {
            CheckDbConsistencyFatal();
            v_helper_ = std::shared_ptr<VHelper>(new VHelper(*this));
        }

        VJHits Align(core::Read& read);

    private:
        DECL_LOGGER("VJQueryAligner");

    public:
        //TODO: make private VHelper JHelper
        class VHelper {
            private:
                VJQueryAligner &vj_query_aligner_;
                CustomGermlineDbHelper kmer_index_helper_;
                algorithms::SubjectQueryKmerIndex<germline_utils::CustomGeneDatabase, seqan::Dna5String> kmer_index_;
                std::shared_ptr<algorithms::PairwiseBlockAligner<germline_utils::CustomGeneDatabase, seqan::Dna5String> > aligner_;

            public:
                VHelper(VJQueryAligner &vj_query_aligner_) :
                    vj_query_aligner_(vj_query_aligner_),
                    kmer_index_helper_(vj_query_aligner_.v_custom_db_),
                    kmer_index_(vj_query_aligner_.v_custom_db_,
                               vj_query_aligner_.algorithm_params_.aligner_params.word_size_v,
                               kmer_index_helper_),
                    aligner_(vj_query_aligner_.get_aligner(
                                kmer_index_, kmer_index_helper_,
                                vj_query_aligner_.CreateBlockAlignmentScoring<VJFinderConfig::AlgorithmParams::ScoringParams::VScoringParams>(
                                    vj_query_aligner_.algorithm_params_.scoring_params.v_scoring),
                                vj_query_aligner_.CreateVBlockAlignerParams())) {
                    //    TRACE("Kmer index for V gene segment DB was constructed");
                }

                std::shared_ptr<algorithms::PairwiseBlockAligner<germline_utils::CustomGeneDatabase, seqan::Dna5String> > get_aligner() const {
                    return aligner_;
                }
        };

        class JHelper {
            private:
                VJQueryAligner &vj_query_aligner_;
                germline_utils::ImmuneGeneType immune_gene_type_;
                const germline_utils::ImmuneGeneDatabase& j_gene_db_;
                ImmuneGeneGermlineDbHelper kmer_index_helper_;
                algorithms::SubjectQueryKmerIndex<germline_utils::ImmuneGeneDatabase, seqan::Dna5String> kmer_index_;
                std::shared_ptr<algorithms::PairwiseBlockAligner<germline_utils::ImmuneGeneDatabase, seqan::Dna5String> > aligner_;

            public:
                JHelper(VJQueryAligner &vj_query_aligner_, germline_utils::ImmuneGeneType immune_gene_type) :
                    vj_query_aligner_(vj_query_aligner_),
                    immune_gene_type_(immune_gene_type),
                    j_gene_db_(vj_query_aligner_.j_custom_db_.GetConstDbByGeneType(immune_gene_type)),
                    kmer_index_helper_(j_gene_db_),
                    kmer_index_(j_gene_db_,
                               vj_query_aligner_.algorithm_params_.aligner_params.word_size_j,
                               kmer_index_helper_),
                    aligner_(vj_query_aligner_.get_aligner(
                                                         kmer_index_, kmer_index_helper_,
                                                         vj_query_aligner_.CreateBlockAlignmentScoring<VJFinderConfig::AlgorithmParams::ScoringParams::JScoringParams>(
                                                                vj_query_aligner_.algorithm_params_.scoring_params.j_scoring),
                                                         vj_query_aligner_.CreateJBlockAlignerParams())) {
                    //    TRACE("Kmer index for J gene segment DB was constructed");
                }

                const germline_utils::ImmuneGeneType& get_immune_gene_type() const {
                    return immune_gene_type_;
                }

                const germline_utils::ImmuneGeneDatabase& get_j_gene_db() const {
                    return j_gene_db_;
                }

                std::shared_ptr<algorithms::PairwiseBlockAligner<germline_utils::ImmuneGeneDatabase, seqan::Dna5String> > get_aligner() {
                    return aligner_;
                }
        };

    private:
        std::shared_ptr<VHelper> v_helper_;
        std::vector<std::shared_ptr<JHelper> > j_helpers_;

        std::shared_ptr<JHelper> get_j_helper(germline_utils::ImmuneGeneType &immune_gene_type) {
            auto it = std::find_if(j_helpers_.begin(), j_helpers_.end(), [&immune_gene_type] (const std::shared_ptr<JHelper> j_helper) {
                    return j_helper->get_immune_gene_type() == immune_gene_type;
                });
            if (it != j_helpers_.end()) {
                return *it;
            }
            j_helpers_.push_back(std::shared_ptr<JHelper>(new JHelper(*this, immune_gene_type)));
            return j_helpers_.back();
        }
    };
}
