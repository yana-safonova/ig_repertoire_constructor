//
// Created by Andrew Bzikadze on 14/05/16.
//

#pragma once

#include "io/library.hpp"
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
        enum class AlignmentCheckerMethod { NoGapsAlignmentCheckerConfig };
        AlignmentCheckerMethod alignment_checker_method;
    };

    struct alignment_cropper_params {
        struct alignment_cropper_method_params { };
        struct upto_reliable_kmer_cropper_params :
            public alignment_cropper_method_params {
            unsigned int kmer_len;
        };
        enum class AlignmentCropperMethod { UptoLastReliableKMer };
        AlignmentCropperMethod alignment_cropper_method;

        upto_reliable_kmer_cropper_params rkmp;
    };

    struct mismatch_finder_params {
        struct mismatch_finder_method_params { };
        struct trivial_mismatch_finder_params:
            public mismatch_finder_method_params { };
        struct no_kneighbours_mismatch_finder_params:
            public mismatch_finder_method_params {
            unsigned int kmer_len;
        };

        enum class MismatchFinderMethod { Trivial, NoKNeighbours };
        MismatchFinderMethod mismatch_finder_method;

        trivial_mismatch_finder_params tmfp;
        no_kneighbours_mismatch_finder_params nknmfp;
    };

    io_params io;
    alignment_checker_params achp;
    alignment_cropper_params acrp;
    mismatch_finder_params mfp;
};

void load(shm_config &cfg, std::string const &filename);

typedef config_common::config<shm_config> shm_cfg;