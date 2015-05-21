#include "ig_config.hpp"
#include "../include/config_common.hpp"

std::string update_hgraph_dir_name(ig_config &cfg) {
	std::stringstream ss;
	ss << path::append_path(cfg.io.output_dir, cfg.hgc_params.hgc_io_params.hg_output_dir);
	ss << "_tau_" << cfg.aligning_params.overlap_mismatches_threshold;
	return ss.str();
}

void AddOutputDir(ig_config::io_params &io) {
	io.load_from = path::append_path(io.output_dir, io.load_from);
	io.log_filename = path::append_path(io.output_dir, io.log_filename);
	io.output_saves = path::append_path(io.output_dir, io.output_saves);
	io.temp_files = path::append_path(io.output_dir, io.temp_files);
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

void load(ig_config::hg_clusterization_params::hg_clusterization_io_params &io_params, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
	load(io_params.hg_output_dir, pt, "hg_output_dir");
    load(io_params.path_to_metis, pt, "path_to_metis");
    load(io_params.run_metis, pt, "run_metis");
    load(io_params.trash_output, pt, "trash_output");
    load(io_params.output_dense_subgraphs, pt, "output_dense_subgraphs");
    load(io_params.dense_subgraphs_dir, pt, "dense_subgraphs_dir");
    io_params.run_metis = path::append_path(io_params.path_to_metis, io_params.run_metis);
}

void load(ig_config::hg_clusterization_params &params, boost::property_tree::ptree const &pt, bool complete) {
    using config_common::load;
    load(params.edge_perc_threshold, pt, "edge_perc_threshold");
    load(params.class_joining_edge_threshold, pt, "class_joining_edge_threshold");
    load(params.min_recessive_abs_size, pt, "min_recessive_abs_size");
    load(params.min_recessive_rel_size, pt, "min_recessive_rel_size");
    load(params.hgc_io_params, pt, "hgc_io_params", complete);
}

void load(ig_config &cfg, boost::property_tree::ptree const &pt, bool complete) {
    using config_common::load;
    load(cfg.io, pt, "io", complete);
    load(cfg.rp, pt, "rp", complete);
    load(cfg.sg, pt, "sg", complete);
    load(cfg.aligning_params, pt, "aligning_params", complete);
    load(cfg.hgc_params, pt, "hgc_params", complete);

    // temporary
    cfg.hgc_params.hgc_io_params.hg_output_dir = update_hgraph_dir_name(cfg);
    cfg.hgc_params.hgc_io_params.dense_subgraphs_dir = path::append_path(cfg.io.output_dir, cfg.hgc_params.hgc_io_params.dense_subgraphs_dir);
    cfg.hgc_params.hgc_io_params.trash_output = path::append_path(cfg.io.temp_files, cfg.hgc_params.hgc_io_params.trash_output);
}

void load(ig_config &cfg, std::string const &filename) {
    boost::property_tree::ptree pt;
    boost::property_tree::read_info(filename, pt);
    load(cfg, pt, true);
}
