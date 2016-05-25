#include "cdr_config.hpp"
#include <config_common.hpp>

namespace cdr_labeler {
    void CDRLabelerConfig::load(std::string config_fname) {
        boost::property_tree::ptree pt;
        boost::property_tree::read_info(config_fname, pt);
        input_params.load(pt, true);
        output_params.load(pt, true);
        run_params.load(pt, true);
        vj_finder::load(vj_finder_config, input_params.vj_finder_config);
    }

    void CDRLabelerConfig::InputParams::load(boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(input_reads, pt, "input_reads");
        load(vj_finder_config, pt, "vj_finder_config");
    }

    void CDRLabelerConfig::OutputParams::load(boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(output_dir, pt, "output_dir");
    }

    void CDRLabelerConfig::RunParams::load(boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(num_threads, pt, "num_threads");
    }
}