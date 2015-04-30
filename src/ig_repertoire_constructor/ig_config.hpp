#pragma once

#include "io/library.hpp"
#include "config_singl.hpp"
#include <boost/property_tree/ptree_fwd.hpp>

struct ig_config {
    struct io_params {
        std::string     log_filename;
        std::string     output_dir;
        io::DataSet<> 	dataset;
        std::string 	output_saves;
        std::string 	load_from;
        std::string     temp_files;
        bool            output_ham_graphs;
        std::string     hgraph_dir;
    };

    struct run_params {
    	bool 			developer_mode;
    	std::string 	entry_point;
    	unsigned        threads_count;
        unsigned        max_memory;
    };

    struct read_aligning_params {
        unsigned min_overlap_length;
        unsigned threshold_for_single_shift;
        unsigned overlap_mismatches_threshold;
    };

    struct hg_clusterization_params {
        double edge_perc_threshold;
        double class_joining_edge_threshold;
    };

    struct singleton_gluer_params {
        size_t k;
        size_t max_kmer_occurences;
        size_t max_distance;
        size_t min_overlap;
        int common_mismatches_threshold;
        int diverse_mismatches_threshold;
    };

    io_params io;
    run_params rp;
    read_aligning_params aligning_params;
    hg_clusterization_params hgc_params;
    singleton_gluer_params sg;
};

void load(ig_config &cfg, std::string const &filename);

typedef config_common::config<ig_config> ig_cfg;
