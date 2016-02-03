#pragma once

#include "io/library.hpp"
#include "config_singl.hpp"
#include <boost/property_tree/ptree_fwd.hpp>
#include "include_me.hpp"

struct dsf_config {
    struct io_params {
        struct input_params {
            std::string 	graph_filename;
        };

        struct output_params {
            std::string     log_filename;
            std::string     output_dir;
            std::string     decomposition_filename;
        };

        struct output_nonparallel_params {
            std::string     graph_copy_filename;
            std::string     permutation_filename;
        };

        struct output_mthreading_params {
            std::string     connected_components_dir;
            std::string     decompositions_dir;
            bool            output_component_decompositions;
        };

        input_params input;
        output_params output_base;
        output_nonparallel_params output_nonparallel;
        output_mthreading_params output_mthreading;
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
        size_t          min_graph_size;
        double          primary_edge_fillin;
        size_t            min_supernode_size;
        double          min_fillin_threshold;
        bool            create_trivial_decomposition;
    };

    io_params io;
    run_params rp;
    dense_sgraph_finder_params dsf_params;
    metis_io_params metis_io;
};

void load(dsf_config &cfg, std::string const &filename);

typedef config_common::config<dsf_config> dsf_cfg;
