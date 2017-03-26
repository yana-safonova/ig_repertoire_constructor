//
// Created by Andrew Bzikadze on 3/15/17.
//

#pragma once

#include "io/library.hpp"
#include <logger/logger.hpp>
#include "config_singl.hpp"
#include "germline_utils/germline_config.hpp"
#include "cdr_config.hpp"

namespace ig_simulator {

struct IgSimulatorConfig {
    struct IOParams {
        struct InputParams {
            germline_utils::GermlineInput germline_input;
            std::string cdr_labeler_config_filename;
        };

        struct OutputParams {
            std::string output_dir;
            std::string log_filename;
        };

        InputParams input_params;
        OutputParams output_params;
    };

    struct AlgorithmParams {
        germline_utils::GermlineParams germline_params;
    };

    IOParams io_params;
    AlgorithmParams algorithm_params;

    cdr_labeler::CDRLabelerConfig cdr_labeler_config;
};

void load(IgSimulatorConfig &cfg, std::string const &filename);

typedef config_common::config<IgSimulatorConfig> igs_cfg;

} // End namespace ig_simulator
