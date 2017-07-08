//
// Created by Andrew Bzikadze on 14/05/16.
//

#pragma once

#include <vector>

#include "config_singl.hpp"
#include <boost/property_tree/ptree_fwd.hpp>

namespace shm_kmer_matrix_estimator {

struct shm_kmer_matrix_estimator_config {
    struct io_params {
        struct input_params {
            std::string v_alignments;
            std::string cdr_details;
        };

        struct output_params {
            std::string output_dir;
            std::string log_filename;
            std::string output_filename_fr;
            std::string output_filename_cdr;
        };

        input_params input;
        output_params output;
    };

    struct alignment_checker_params {
        enum class AlignmentCheckerMethod { NoGaps };
        const static std::vector<std::string> alignment_checker_method_names;
        AlignmentCheckerMethod alignment_checker_method;

        enum class FunctionalityMethod { all, productive, nonproductive };
        const static std::vector<std::string> functionality_methods_names;
        FunctionalityMethod functionality_method;
    };

    struct alignment_cropper_params {
        struct alignment_cropper_method_params {};
        struct upto_reliable_kmer_cropper_params:
            public alignment_cropper_method_params {
            unsigned int kmer_len;
            unsigned int hash_base;
        };
        enum class AlignmentCropperMethod { UptoLastReliableKMer };
        const static std::vector<std::string> alignment_cropper_method_names;
        AlignmentCropperMethod alignment_cropper_method;

        upto_reliable_kmer_cropper_params rkmp;
    };

    struct mutations_strategy_params {
        struct mutations_strategy_method_params {};
        struct trivial_mutations_strategy_params:
            public mutations_strategy_method_params {
        };
        struct no_kneighbours_mutations_strategy_params:
            public mutations_strategy_method_params {
        };

        enum class MutationsStrategyMethod { Trivial, NoKNeighbours };
        const static std::vector<std::string> mutation_strategy_method_names;
        MutationsStrategyMethod mutations_strategy_method;

        trivial_mutations_strategy_params tmfp;
        no_kneighbours_mutations_strategy_params nknmfp;
        unsigned int kmer_len;
    };

    io_params io;
    alignment_checker_params achp;
    alignment_cropper_params acrp;
    mutations_strategy_params mfp;
};

std::istream &operator>>(std::istream &in,
                         shm_kmer_matrix_estimator_config::mutations_strategy_params::MutationsStrategyMethod &strategy);

void load(shm_kmer_matrix_estimator_config &cfg, std::string const &filename);

typedef config_common::config<shm_kmer_matrix_estimator_config> shm_cfg;

using FunctionalityMethod = shm_kmer_matrix_estimator_config::alignment_checker_params::FunctionalityMethod;

} // End namespace shm_kmer_matrix_estimator
