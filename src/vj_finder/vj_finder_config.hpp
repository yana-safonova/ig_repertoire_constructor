#pragma once

#include "io/library.hpp"
#include "config_singl.hpp"
#include <boost/property_tree/ptree_fwd.hpp>
#include "germline_utils/germline_config.hpp"

namespace vj_finder {
    struct VJFinderConfig {
        struct RunParams {
            size_t num_threads;
        };

        struct IOParams {
            struct InputParams {
                std::string input_reads;
                std::string config_dir;
                germline_utils::GermlineInput germline_input;
            };

            struct OutputParams {
                struct OutputFiles {
                    std::string output_dir;
                    std::string log_filename;
                    std::string cleaned_reads_fname;
                    std::string filtered_reads_fname;
                    std::string alignment_info_fname;
                    std::string filtering_info_filename;
                    std::string valignments_filename;
                };

                struct OutputDetails {
                    bool compress;
                    bool fix_spaces;
                    std::string separator;
                    size_t num_aligned_candidates;
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
                bool enable_filtering;
                size_t min_v_segment_length;
                size_t min_j_segment_length;
                int left_uncovered_limit;
                int right_uncovered_limit; // It should be at least 4 (=1 + 3cropped) 1bp trimming is common
                size_t min_aligned_length;
            };

            struct FixCropFillParams {
                enum FixCropFillAlgorithm { UnknowmFCFAlgorithm, AggressiveFCFAlgorithm };

                bool enable_fixing;
                size_t fix_left;
                size_t fix_right;
                bool enable_filling;
                bool fill_left;
                bool fill_right;
                bool enable_cropping;
                bool crop_left;
                bool crop_right;
                FixCropFillAlgorithm fcf_algorithm;
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
            germline_utils::GermlineParams germline_params;
            FilteringParams filtering_params;
            FixCropFillParams fix_crop_fill_params;
            ScoringParams scoring_params;
        };

        RunParams run_params;
        IOParams io_params;
        AlgorithmParams algorithm_params;
    };

    void load(VJFinderConfig &cfg, std::string const &filename);

    typedef config_common::config<VJFinderConfig> vjf_cfg;
    
    void update_output_files_config(VJFinderConfig::IOParams::OutputParams::OutputFiles & of);
}
