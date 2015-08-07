#pragma once

#include "io/library.hpp"
#include "config_singl.hpp"
#include <boost/property_tree/ptree_fwd.hpp>

struct dsf_config {
    struct io_params {
        std::string     log_filename;
        std::string     output_dir;
        std::string 	graph_filename;
        std::string     temp_dir;
    };

    struct run_params {
    	bool 			developer_mode;
    	unsigned        threads_count;
        unsigned        max_memory;
    };

    struct dense_sgraph_finder_params {
        double edge_perc_threshold;
        double class_joining_edge_threshold;

        struct hg_clusterization_io_params {
        	std::string 	hg_output_dir;
            std::string 	path_to_metis;
            std::string		run_metis;
            std::string		trash_output;
            bool			output_dense_subgraphs;
            std::string 	dense_subgraphs_dir;
        };

        hg_clusterization_io_params hgc_io_params;
    };

    io_params io;
    run_params rp;
    dense_sgraph_finder_params hgc_params;
};

void load(dsf_config &cfg, std::string const &filename);

typedef config_common::config<dsf_config> dsf_cfg;
