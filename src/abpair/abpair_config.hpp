#pragma once

#include "io/library.hpp"
#include "config_singl.hpp"
#include <boost/property_tree/ptree_fwd.hpp>

struct abpair_config {
    struct io_config {
        struct input_config {
            std::string input_sequences;
        };

        struct output_config {
            std::string output_dir;
            std::string log_filename;
            bool output_statistics;
            std::string hc_ambiguous_dir;
            bool output_barcodes;
            std::string barcode_dir;
        };

        input_config input;
        output_config output;
    };

    io_config io;
};

void load(abpair_config &cfg, const std::string &filename);

typedef config_common::config<abpair_config> abp_cfg;
