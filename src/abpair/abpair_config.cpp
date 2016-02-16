#include "abpair_config.hpp"
#include "../include/config_common.hpp"

void load(abpair_config::io_config::input_config &input, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(input.input_sequences, pt, "input_sequences");
}

void load(abpair_config::io_config::output_config &output, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(output.output_dir, pt, "output_dir");
    load(output.log_filename, pt, "log_filename");
    load(output.output_statistics, pt, "output_statistics");
    load(output.hc_ambiguous_dir, pt, "hc_ambiguous_dir");
    output.hc_ambiguous_dir = path::append_path(output.output_dir, output.hc_ambiguous_dir);
    load(output.output_barcodes, pt, "output_barcodes");
    load(output.barcode_dir, pt, "barcode_dir");
    output.barcode_dir = path::append_path(output.output_dir, output.barcode_dir);
}

void load(abpair_config::io_config &io, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(io.input, pt, "input");
    load(io.output, pt, "output");
}

void load(abpair_config &cfg, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(cfg.io, pt, "io");
}

void load(abpair_config &cfg, const std::string &filename) {
    boost::property_tree::ptree pt;
    boost::property_tree::read_info(filename, pt);
    load(cfg, pt, true);
}