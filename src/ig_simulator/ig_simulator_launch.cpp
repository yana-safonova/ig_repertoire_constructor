//
// Created by Andrew Bzikadze on 3/15/17.
//

#include "ig_simulator_launch.hpp"
#include "germline_db_generator.hpp"

using namespace vj_finder;

namespace ig_simulator {

void IgSimulatorLaunch::Run() {
    INFO("== IgSimulator starts ==");

    GermlineDbGenerator db_generator(config_.io_params.input_params.germline_input,
                                     config_.algorithm_params.germline_params);
    INFO("Generation of DB for variable segments...");
    germline_utils::CustomGeneDatabase v_db = db_generator.GenerateVariableDb();
    INFO("Generation of DB for join segments...");
    germline_utils::CustomGeneDatabase j_db = db_generator.GenerateJoinDb();
}

} // End namespace ig_simulator
