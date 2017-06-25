#pragma once

#include "vj_finder_config.hpp"
#include <germline_utils/germline_databases/custom_gene_database.hpp>
#include "vj_alignment_structs.hpp"
#include "block_alignment/pairwise_block_aligner.hpp"

using namespace algorithms;
namespace vj_finder {

    class VJQueryAligner {
        struct VHelper {
            VJQueryAligner &vj_query_aligner;
            CustomGermlineDbHelper kmer_index_helper;
            SubjectQueryKmerIndex<germline_utils::CustomGeneDatabase, seqan::Dna5String> kmer_index;
            std::shared_ptr<PairwiseBlockAligner<germline_utils::CustomGeneDatabase, seqan::Dna5String> > aligner;
            

            VHelper(VJQueryAligner &vj_query_aligner_) :
                                vj_query_aligner(vj_query_aligner_),
                                kmer_index_helper(vj_query_aligner.v_custom_db_),
                                kmer_index(vj_query_aligner.v_custom_db_, 
                                           vj_query_aligner.algorithm_params_.aligner_params.word_size_v,
                                           kmer_index_helper),
                                aligner(vj_query_aligner.get_aligner(
                                            kmer_index, kmer_index_helper,
                                            vj_query_aligner.CreateBlockAlignmentScoring<VJFinderConfig::AlgorithmParams::ScoringParams::VScoringParams>(vj_query_aligner.algorithm_params_.scoring_params.v_scoring),
                                            vj_query_aligner.CreateVBlockAlignerParams())) { 
                                //    TRACE("Kmer index for V gene segment DB was constructed");
                                }
        };

        std::shared_ptr<VHelper> v_helper;
        const VJFinderConfig::AlgorithmParams & algorithm_params_;

        struct JHelper {
            VJQueryAligner &vj_query_aligner;
            const germline_utils::ImmuneGeneDatabase& j_gene_db;
            germline_utils::ImmuneGeneType immune_gene_type;
            ImmuneGeneGermlineDbHelper kmer_index_helper;
            SubjectQueryKmerIndex<germline_utils::ImmuneGeneDatabase, seqan::Dna5String> kmer_index;
            std::shared_ptr<PairwiseBlockAligner<germline_utils::ImmuneGeneDatabase, seqan::Dna5String> > aligner; 
        
            

            JHelper(VJQueryAligner &vj_query_aligner_, germline_utils::ImmuneGeneType immune_gene_type) :
                                vj_query_aligner(vj_query_aligner_),
                                immune_gene_type(immune_gene_type),
                                j_gene_db(vj_query_aligner.j_custom_db_.GetConstDbByGeneType(immune_gene_type)),
                                kmer_index_helper(j_gene_db),
                                kmer_index(j_gene_db, 
                                           vj_query_aligner.algorithm_params_.aligner_params.word_size_j, 
                                           kmer_index_helper),
                                aligner(vj_query_aligner.get_aligner(
                                            kmer_index, kmer_index_helper,
                                            vj_query_aligner.CreateBlockAlignmentScoring<VJFinderConfig::AlgorithmParams::ScoringParams::JScoringParams>(
                                            vj_query_aligner.algorithm_params_.scoring_params.j_scoring),
                                            vj_query_aligner.CreateJBlockAlignerParams())) { 
                                //    TRACE("Kmer index for J gene segment DB was constructed");
                                }
        };

        std::vector<std::shared_ptr<JHelper> > j_helpers;

        std::shared_ptr<JHelper> get_j_helper(germline_utils::ImmuneGeneType &immune_gene_type) {
            auto it = std::find_if(j_helpers.begin(), j_helpers.end(), [&immune_gene_type] (const std::shared_ptr<JHelper> j_helper) {
                return j_helper->immune_gene_type == immune_gene_type;
            });
            if (it != j_helpers.end()) {
                return *it;
            }
            j_helpers.push_back(std::shared_ptr<JHelper>(new JHelper(*this, immune_gene_type)));
            return j_helpers.back();
        }

        core::ReadArchive &read_archive_;
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
            switch(algorithm_params_.aligner_params.aligner_algorithm) {
            case vj_finder::VJFinderConfig::AlgorithmParams::AlignerParams::AlignerAlgorithm::QuadraticDAGAlignerAlgorithm:
                return std::shared_ptr<algorithms::PairwiseBlockAligner<SubjectDatabase, StringType> >(
                    new algorithms::QuadraticDAGPairwiseBlockAligner<SubjectDatabase, StringType>(kmer_index, kmer_index_helper, scoring, params));
            
            case vj_finder::VJFinderConfig::AlgorithmParams::AlignerParams::AlignerAlgorithm::LisAlignerAlgorithm:
                return std::shared_ptr<algorithms::PairwiseBlockAligner<SubjectDatabase, StringType> >(
                    new algorithms::LisPairwiseBlockAligner<SubjectDatabase, StringType>(kmer_index, kmer_index_helper, scoring, params));        
            
            case vj_finder::VJFinderConfig::AlgorithmParams::AlignerParams::AlignerAlgorithm::QuadraticDpAlignerAlgorithm:
                return std::shared_ptr<algorithms::PairwiseBlockAligner<SubjectDatabase, StringType> >(
                    new algorithms::QuadraticDpPairwiseBlockAligner<SubjectDatabase, StringType>(kmer_index, kmer_index_helper, scoring, params));        
            
            default:
                VERIFY_MSG(false, "Failed to determine block alignment algorithm, using the default one");
            };
            return NULL;
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
            v_helper = std::shared_ptr<VHelper>(new VHelper(*this));
        }

        VJHits Align(const core::Read& read);

    private:
        DECL_LOGGER("VJQueryAligner");
    };
}
