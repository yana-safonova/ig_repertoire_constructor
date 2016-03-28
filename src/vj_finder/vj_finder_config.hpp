#pragma once

#include "io/library.hpp"
#include "config_singl.hpp"
#include <boost/property_tree/ptree_fwd.hpp>

namespace vj_finder {
    struct vjf_config {
        struct run_params {
            size_t num_threads;
        };

        struct io_params {
            struct input_params {
                std::string input_reads;
            };

            struct output_params {
                struct output_files {
                    std::string output_dir;
                    std::string output_filename;
                    std::string bad_output_filename;
                    std::string add_info_filename;
                    std::string discard_info_filename;
                    std::string valignments_filename;
                };

                struct output_details {
                    bool verbose;
                    bool compress;
                    bool fix_spaces;
                    std::string separator;
                };

                output_files of;
                output_details od;
            };

            input_params input;
            output_params output;
        };

        struct algorithm_params {
            struct aligner_params {
                size_t left_uncovered_limit;
                size_t right_uncovered_limit; // It should be at least 4 (=1 + 3cropped) 1bp trimming is common
                size_t min_v_segment_length;
                size_t min_j_segment_length;
                int K;
                int word_size_j;
                int min_k_coverage;
                int min_k_coverage_j;

            };

            struct germline_params {
                std::string organism;
                std::string loci;
                std::string db_directory;
                bool pseudogenes;
            };

            struct query_params {
                int max_candidates;
                int max_candidates_j;
                bool fix_strand;
                bool consistent_loci;
                size_t min_len;
            };

            struct fix_crop_fill_params {
                bool fill_left;
                bool fill_right;
                size_t fix_left;
                size_t fix_right;
                bool crop_left;
                bool crop_right;
            };

            struct alignment_scoring_params {
                int max_global_gap = 24;
                int max_local_deletions = 12;
                int max_local_insertions = 12;

                int gap_opening_cost = 4;
                int gap_extention_cost = 1;
            };

            aligner_params ap;
            germline_params gp;
            query_params qp;
            fix_crop_fill_params fxp;
            alignment_scoring_params asp;
        };

        run_params rp;
        io_params iop;
        algorithm_params algop;
    };

    void load(vjf_config &cfg, std::string const &filename);
    typedef config_common::config<vjf_config> vjf_cfg;
}
