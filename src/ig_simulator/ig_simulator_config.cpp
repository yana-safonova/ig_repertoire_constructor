//
// Created by Andrew Bzikadze on 3/15/17.
//

#include "ig_simulator_config.hpp"
#include <boost/property_tree/ptree_fwd.hpp>
#include <config_common.hpp>

namespace ig_simulator {

void load(IgSimulatorConfig::IOParams::InputParams &input_params, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(input_params.germline_input, pt, "germline_input");
    load(input_params.cdr_labeler_config_filename, pt, "cdr_labeler_config_filename");
}

void load(IgSimulatorConfig::IOParams::OutputParams &output_params, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(output_params.log_filename, pt, "log_filename");
    load(output_params.output_dir, pt, "output_dir");
}

void load(IgSimulatorConfig::IOParams &io_params, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(io_params.input_params, pt, "input_params");
    load(io_params.output_params, pt, "output_params");
}

void load(IgSimulatorConfig::AlgorithmParams &algorithm_params,
          boost::property_tree::ptree const &pt, bool)
{
    using config_common::load;
    load(algorithm_params.germline_params, pt, "germline_params");
}

void load(IgSimulatorConfig &cfg, boost::property_tree::ptree const &pt, bool complete) {
    using config_common::load;
    load(cfg.io_params, pt, "io_params", complete);
    load(cfg.algorithm_params, pt, "algorithm_params", complete);
    cfg.cdr_labeler_config.load(cfg.io_params.input_params.cdr_labeler_config_filename);
}

void load(IgSimulatorConfig &cfg, std::string const &filename) {
    boost::property_tree::ptree pt;
    boost::property_tree::read_info(filename, pt);
    load(cfg, pt, true);
}

} // End namespace ig_simulator
