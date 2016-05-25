#pragma once

#include "vj_finder_config.hpp"

namespace cdr_labeler {
    struct CDRLabelerConfig {
        struct InputParams {
            std::string input_reads;
            std::string vj_finder_config;

            void load(boost::property_tree::ptree const &pt, bool complete);
        };

        struct OutputParams {
            std::string output_dir;

            void load(boost::property_tree::ptree const &pt, bool complete);
        };

        struct RunParams {
            size_t num_threads;

            void load(boost::property_tree::ptree const &pt, bool complete);
        };

        InputParams input_params;
        OutputParams output_params;
        RunParams run_params;
        vj_finder::VJFinderConfig vj_finder_config;

        void load(std::string config_fname);
    };
}
