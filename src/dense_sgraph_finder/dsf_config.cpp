#include "dsf_config.hpp"
#include "../include/config_common.hpp"

void updateIO(dsf_config::io_params &io) {
    io.log_filename = path::append_path(io.output_dir, io.log_filename);
    io.decomposition_filename = path::append_path(io.output_dir, io.decomposition_filename);
}

void updateMetisIO(dsf_config::io_params &io, dsf_config::metis_io_params &metis_io) {
    metis_io.trash_output = path::append_path(io.output_dir, metis_io.trash_output);
}

void load(dsf_config::io_params &io, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(io.log_filename, pt, "log_filename");
    load(io.output_dir, pt, "output_dir");
    load(io.graph_filename, pt, "graph_filename");
    load(io.decomposition_filename, pt, "decomposition_filename");
    updateIO(io);
}

void load(dsf_config::run_params &rp, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(rp.developer_mode, pt, "developer_mode");
    load(rp.threads_count, pt, "threads_count");
    load(rp.max_memory, pt, "max_memory");
}

void load(dsf_config::dense_sgraph_finder_params &params, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(params.edge_perc_threshold, pt, "edge_perc_threshold");
    load(params.class_joining_edge_threshold, pt, "class_joining_edge_threshold");
}

void load(dsf_config::metis_io_params &metis_io, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(metis_io.path_to_metis, pt, "path_to_metis");
    load(metis_io.run_metis, pt, "run_metis");
    load(metis_io.trash_output, pt, "trash_output");
    metis_io.run_metis = path::append_path(metis_io.path_to_metis, metis_io.run_metis);
}

void load(dsf_config &cfg, boost::property_tree::ptree const &pt, bool complete) {
    using config_common::load;
    load(cfg.io, pt, "io", complete);
    load(cfg.rp, pt, "rp", complete);
    load(cfg.dsf_params, pt, "dsf_params", complete);
    load(cfg.metis_io, pt, "metis_io", complete);
    updateMetisIO(cfg.io, cfg.metis_io);
}

void load(dsf_config &cfg, std::string const &filename) {
    boost::property_tree::ptree pt;
    boost::property_tree::read_info(filename, pt);
    load(cfg, pt, true);
}
