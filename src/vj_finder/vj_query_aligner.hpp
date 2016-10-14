#pragma once

#include "vj_finder_config.hpp"
#include <germline_utils/germline_databases/custom_gene_database.hpp>
#include "vj_alignment_structs.hpp"
#include "block_alignment/pairwise_block_aligner.hpp"

namespace vj_finder {
    class VJQueryAligner {
        const VJFinderConfig::AlgorithmParams & algorithm_params_;

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