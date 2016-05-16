#pragma once

#include "io/library.hpp"
#include "config_singl.hpp"
#include <boost/property_tree/ptree_fwd.hpp>

namespace vj_finder {
    struct vjf_config {
        struct RunParams {
            size_t num_threads;
        };

        struct IOParams {
            struct InputParams {
                struct GermlineInput {
                    std::string ig_dir;
                    std::string tcr_dir;
                    std::string germline_filenames_config;
                };

                std::string input_reads;
                std::string config_dir;
                GermlineInput germline_input;
            };

            struct OutputParams {
                struct OutputFiles {
                    std::string output_dir;
                    std::string log_filename;
                    std::string output_filename;
                    std::string bad_output_filename;
                    std::string add_info_filename;
                    std::string discard_info_filename;
                    std::string valignments_filename;
                };

                struct OutputDetails {
                    bool verbose;
                    bool compress;
                    bool fix_spaces;
                    std::string separator;
                };

                OutputFiles output_files;
                OutputDetails output_details;
            };

            InputParams input_params;
            OutputParams output_params;
        };

        struct AlgorithmParams {
            struct AlignerParams {
                size_t word_size_v;
                size_t word_size_j;
                size_t min_k_coverage_v;
                size_t min_k_coverage_j;
                size_t max_candidates_v;
                size_t max_candidates_j;
                bool fix_strand;
            };

            struct FilteringParams {
                size_t min_v_segment_length;
                size_t min_j_segment_length;
                size_t left_uncovered_limit;
                size_t right_uncovered_limit; // It should be at least 4 (=1 + 3cropped) 1bp trimming is common
                size_t min_aligned_length;
            };

            struct GermlineParams {
                std::string germline_dir;
                std::string organism;
                std::string loci;
                bool pseudogenes;
            };

            struct FixCropFillParams {
                bool fill_left;
                bool fill_right;
                size_t fix_left;
                size_t fix_right;
                bool crop_left;
                bool crop_right;
            };

            struct ScoringParams {
                struct VScoringParams {
                    int max_global_gap;
                    int max_local_insertions;
                    int max_local_deletions;
                    int gap_opening_cost;
                    int gap_extention_cost;
                    int mismatch_extention_cost;
                    int mismatch_opening_cost;
                    int match_reward;
                };

                struct JScoringParams {
                    int max_global_gap;
                    int max_local_insertions;
                    int max_local_deletions;
                    int gap_opening_cost;
                    int gap_extention_cost;
                    int mismatch_extention_cost;
                    int mismatch_opening_cost;
                    int match_reward;
                };

                VScoringParams v_scoring;
                JScoringParams j_scoring;
            };

            AlignerParams aligner_params;
            GermlineParams germline_params;
            FilteringParams filtering_params;
            FixCropFillParams fix_crop_fill_params;
            ScoringParams scoring_params;
        };

        RunParams run_params;
        IOParams io_params;
        AlgorithmParams algorithm_params;
    };

    void load(vjf_config &cfg, std::string const &filename);

    typedef config_common::config<vjf_config> vjf_cfg;
}
