#pragma once

#include <cdr_config.hpp>

namespace antevolo {
    struct AntEvoloConfig {
        struct InputParams {
            std::string input_reads;
            std::string cdr_labeler_config_fname;

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
            std::string tree_shms;
        };

        struct RunParams {
            int num_threads;
        };

        struct AlgorithmParams {
            struct SimilarCDR3Params {
                size_t num_mismatches;
                size_t num_indels;
            };

            struct EdgeConstructionParams {
                double intersected_edge_coeff;
                size_t min_num_intersected_v_shms;
            };

            SimilarCDR3Params similar_cdr3s_params;
            EdgeConstructionParams edge_construction_params;
        };

        InputParams input_params;
        OutputParams output_params;
        RunParams run_params;
        AlgorithmParams algorithm_params;
        cdr_labeler::CDRLabelerConfig cdr_labeler_config;

        void load(std::string antevolo_config_fname);
    };
}