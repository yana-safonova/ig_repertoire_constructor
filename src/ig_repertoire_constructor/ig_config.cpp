#include "ig_config.hpp"
#include "../include/config_common.hpp"

std::string update_hgraph_dir_name(std::string hgraph_dir, size_t overlap_mismatches_threshold) {
    std::stringstream ss;
    ss << hgraph_dir << "_" << overlap_mismatches_threshold;
    return ss.str();
}

void load(ig_config::io_params &io, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(io.log_filename, pt, "log_filename");
    load(io.output_dir, pt, "output_dir");
    load(io.output_saves, pt, "output_saves");
    load(io.load_from, pt, "load_from");
    load(io.temp_files, pt, "temp_files");

    std::string fname;
    load(fname, pt, "dataset");
    io.dataset.load(fname);

    load(io.output_ham_graphs, pt, "output_ham_graphs");
    load(io.hgraph_dir, pt, "hgraph_dir");
    io.hgraph_dir = io.output_saves + "/" + io.hgraph_dir;
}

void load(ig_config::run_params &rp, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(rp.developer_mode, pt, "developer_mode");
    load(rp.entry_point, pt, "entry_point");
    load(rp.threads_count, pt, "threads_count");
    load(rp.max_memory, pt, "max_memory");
}

void load(ig_config::singleton_gluer_params &sg, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(sg.k, pt, "k");
    load(sg.max_kmer_occurences, pt, "max_kmer_occurences");
    load(sg.max_distance, pt, "max_distance");
    load(sg.min_overlap, pt, "min_overlap");
    load(sg.common_mismatches_threshold, pt, "common_mismatches_threshold");
    load(sg.diverse_mismatches_threshold, pt, "diverse_mismatches_threshold");
}

void load(ig_config::read_aligning_params &params, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(params.min_overlap_length, pt, "min_overlap_length");
    load(params.threshold_for_single_shift, pt, "threshold_for_single_shift");
    load(params.overlap_mismatches_threshold, pt, "overlap_mismatches_threshold");
}

void load(ig_config::hg_clusterization_params &params, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(params.edge_perc_threshold, pt, "edge_perc_threshold");
    load(params.class_joining_edge_threshold, pt, "class_joining_edge_threshold");
}

void load(ig_config &cfg, boost::property_tree::ptree const &pt, bool complete) {
    using config_common::load;
    load(cfg.io, pt, "io", complete);
    load(cfg.rp, pt, "rp", complete);
    load(cfg.sg, pt, "sg", complete);
    load(cfg.aligning_params, pt, "aligning_params", complete);
    load(cfg.hgc_params, pt, "hgc_params");

    // temporary
    cfg.io.hgraph_dir = update_hgraph_dir_name(cfg.io.hgraph_dir, cfg.aligning_params.overlap_mismatches_threshold);
}

void load(ig_config &cfg, std::string const &filename) {
    boost::property_tree::ptree pt;
    boost::property_tree::read_info(filename, pt);
    load(cfg, pt, true);
}
