#pragma once

#include "vj_finder_config.hpp"

namespace cdr_labeler {
    struct CDRLabelerConfig {
        struct InputParams {
            std::string input_reads;
            std::string vj_finder_config;
        };

        struct OutputParams {
            std::string output_dir;
        };

        struct RunParams {
            size_t num_threads;
        };

        struct CDRsParams {
            struct HCDR1Params {
                size_t start_pos;
                size_t start_shift;
                size_t min_length;
                size_t max_length;
                std::string residues_before;
                std::string residues_after;
            };

            struct HCDR2Params {
                size_t distance_from_cdr1_end;
                size_t distance_shift;
                size_t min_length;
                size_t max_length;
                std::string residues_before;
                std::string residues_after;
            };

            struct HCDR3Params {
                size_t distance_from_cdr2_end;
                size_t distance_shift;
                size_t min_length;
                size_t max_length;
                std::string residues_before;
                std::string residues_after;
            };

            HCDR1Params hcdr1_params;
            HCDR2Params hcdr2_params;
            HCDR3Params hcdr3_params;
        };

        InputParams input_params;
        OutputParams output_params;
        RunParams run_params;
        CDRsParams cdrs_params;
        vj_finder::VJFinderConfig vj_finder_config;

        void load(std::string config_fname);
    };
}
