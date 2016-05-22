//
// Created by Andrew Bzikadze on 14/05/16.
//

#pragma once

#include "config_singl.hpp"
#include <boost/property_tree/ptree_fwd.hpp>

struct shm_config {
    struct io_params {
        struct input_params {
            std::string input_filename;
        };

        struct output_params {
            std::string output_dir;
            std::string log_filename;
            std::string output_filename;
        };

        input_params input;
        output_params output;
    };

    struct alignment_checker_params {
        enum class AlignmentCheckerMethod { NoGapsAlignmentChecker};
        AlignmentCheckerMethod alignment_checker_method;
    };

    struct alignment_cropper_params {
        struct alignment_cropper_method_params { };
        struct upto_reliable_kmer_cropper_params :
            public alignment_cropper_method_params {
            unsigned int kmer_len;
            unsigned int hash_base;
        };
        enum class AlignmentCropperMethod { UptoLastReliableKMer };
        AlignmentCropperMethod alignment_cropper_method;

        upto_reliable_kmer_cropper_params rkmp;
    };

    struct mutations_strategy_params {
        struct mutations_strategy_method_params { };
        struct trivial_mutations_strategy_params:
            public mutations_strategy_method_params { };
        struct no_kneighbours_mutations_strategy_params:
            public mutations_strategy_method_params {
            unsigned int kmer_len;
        };

        enum class MutationsStrategyMethod { Trivial, NoKNeighbours };
        MutationsStrategyMethod mutations_strategy_method;

        trivial_mutations_strategy_params tmfp;
        no_kneighbours_mutations_strategy_params nknmfp;
    };

    io_params io;
    alignment_checker_params achp;
    alignment_cropper_params acrp;
    mutations_strategy_params mfp;
};

void load(shm_config &cfg, std::string const &filename);

typedef config_common::config<shm_config> shm_cfg;