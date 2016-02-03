#include "dsf_config.hpp"
#include "../include/config_common.hpp"

void updateIO(dsf_config::io_params &io) {
    io.output_base.log_filename = path::append_path(io.output_base.output_dir, io.output_base.log_filename);
    io.output_base.decomposition_filename = path::append_path(io.output_base.output_dir,
                                                              io.output_base.decomposition_filename);
    io.output_nonparallel.graph_copy_filename = path::append_path(io.output_base.output_dir,
                                                                  io.output_nonparallel.graph_copy_filename);
    io.output_nonparallel.permutation_filename = path::append_path(io.output_base.output_dir,
                                                                   io.output_nonparallel.permutation_filename);
    io.output_mthreading.connected_components_dir = path::append_path(io.output_base.output_dir,
                                                                      io.output_mthreading.connected_components_dir);
    io.output_mthreading.decompositions_dir = path::append_path(io.output_base.output_dir,
                                                                      io.output_mthreading.decompositions_dir);
}

void updateMetisIO(dsf_config::io_params &io, dsf_config::metis_io_params &metis_io) {
    metis_io.trash_output = path::append_path(io.output_base.output_dir, metis_io.trash_output);
}

void load(dsf_config::io_params::input_params &input_params, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(input_params.graph_filename, pt, "graph_filename");
}

void load(dsf_config::io_params::output_params &output_base, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(output_base.log_filename, pt, "log_filename");
    load(output_base.output_dir, pt, "output_dir");
    load(output_base.decomposition_filename, pt, "decomposition_filename");
}

void load(dsf_config::io_params::output_nonparallel_params &output_nonparallel,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(output_nonparallel.graph_copy_filename, pt, "graph_copy_filename");
    load(output_nonparallel.permutation_filename, pt, "permutation_filename");
}

void load(dsf_config::io_params::output_mthreading_params &output_mthreading,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(output_mthreading.connected_components_dir, pt, "connected_components_dir");
    load(output_mthreading.decompositions_dir, pt, "decompositions_dir");
    load(output_mthreading.output_component_decompositions, pt, "output_component_decompositions");
}

void load(dsf_config::io_params &io, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(io.input, pt, "input");
    load(io.output_base, pt, "output_base");
    load(io.output_nonparallel, pt, "output_nonparallel");
    load(io.output_mthreading, pt, "output_mthreading");
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
    load(params.min_graph_size, pt, "min_graph_size");
    load(params.min_fillin_threshold, pt, "min_fillin_threshold");
    load(params.min_supernode_size, pt, "min_supernode_size");
    load(params.primary_edge_fillin, pt, "primary_edge_fillin");
    load(params.create_trivial_decomposition, pt, "create_trivial_decomposition");
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
