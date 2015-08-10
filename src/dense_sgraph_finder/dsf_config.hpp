#pragma once

#include "io/library.hpp"
#include "config_singl.hpp"
#include <boost/property_tree/ptree_fwd.hpp>
#include "include_me.hpp"

struct dsf_config {
    struct io_params {
        std::string     log_filename;
        std::string     output_dir;
        std::string 	graph_filename;
        std::string     decomposition_filename;
    };

    struct run_params {
    	bool 			developer_mode;
    	unsigned        threads_count;
        unsigned        max_memory;
    };

    struct metis_io_params {
        std::string 	path_to_metis;
        std::string		run_metis;
        std::string		trash_output;
    };

    struct dense_sgraph_finder_params {
        double edge_perc_threshold;
        double class_joining_edge_threshold;
    };

    io_params io;
    run_params rp;
    dense_sgraph_finder_params dsf_params;
    metis_io_params metis_io;
};

void load(dsf_config &cfg, std::string const &filename);

typedef config_common::config<dsf_config> dsf_cfg;
