#pragma once

#include <cdr_config.hpp>
#include <shm_kmer_matrix_estimator_config.hpp>
#include <germline_utils/chain_type.hpp>

namespace antevolo {
    struct AntEvoloConfig {
        struct InputParams {
            std::string input_reads;
            std::string decomposition_rcm;
            std::string cdr_labeler_config_fname;
            std::string shm_kmer_matrix_estimator_config_fname;

            std::string shm_kmer_model_igh;
            std::string shm_kmer_model_igk;
            std::string shm_kmer_model_igl;
        };

        struct OutputParams {
            std::string output_dir;
            std::string cdr_graph_dir;
            std::string tree_dir;
            std::string vertex_dir;
            std::string trash_output;

            std::string tree_details;
            std::string tree_shm_dir;

            struct ParallelSHMOutput {
                std::string parallel_bulges_dir;
                std::string parallel_shm_dir;
            };

            ParallelSHMOutput parallel_shm_output;
        };

        struct RunParams {
            int num_threads;
        };

        struct AlgorithmParams {
            bool compare;
            bool model;

            struct SimilarCDR3Params {
                size_t num_mismatches_igh;
                size_t num_mismatches_igk;
                size_t num_mismatches_igl;
                size_t num_indels;
            };

            struct EdgeConstructionParams {
                double intersected_edge_coeff;
                size_t min_num_intersected_v_shms;
            };

            struct ParallelEvolutionParams {
                bool enable_parallel_shms_finder;
            };

            SimilarCDR3Params similar_cdr3s_params;
            EdgeConstructionParams edge_construction_params;
            ParallelEvolutionParams parallel_evolution_params;

            size_t GetNumMismatchesByChainType(const germline_utils::ImmuneChainType chain) const;
        };

        InputParams input_params;
        OutputParams output_params;
        RunParams run_params;
        AlgorithmParams algorithm_params;
        cdr_labeler::CDRLabelerConfig cdr_labeler_config;
        shm_kmer_matrix_estimator::shm_kmer_matrix_estimator_config shm_config;

        void load(std::string antevolo_config_fname);
    };
}