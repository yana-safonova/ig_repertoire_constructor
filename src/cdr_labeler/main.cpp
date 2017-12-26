#include <logger/logger.hpp>
#include <logger/log_writers.hpp>
#include <verify.hpp>
#include <segfault_handler.hpp>
#include <perfcounter.hpp>

#include <copy_file.hpp>

#include "cdr_config.hpp"
#include "cdr_launch.hpp"
#include "../config/config_utils.hpp"
#include "CdrLConfigLoader.hpp"

namespace {
    void create_console_logger(/*std::string cfg_filename*/) {
        using namespace logging;
        std::string log_props_file = "";
        //if (!path::FileExists(log_props_file)){
        //    log_props_file = path::append_path(path::parent_path(cfg_filename), log_props_file);
        //}
        logger *lg = create_logger(path::FileExists(log_props_file) ? log_props_file : "");
        lg->add_writer(std::make_shared<console_writer>());
        attach_logger(lg);
    }
}

int main(int argc, char **argv) {
    omp_set_num_threads(1);
    segfault_handler sh;
    perf_counter pc;
    create_console_logger();
    cdr_labeler::CdrLConfigLoader().LoadConfig(argc, argv);
    const cdr_labeler::CDRLabelerConfig config = cdr_labeler::cdrl_cfg::get();
    cdr_labeler::CDRLabelerLaunch(config).Launch();
    //create_console_logger(cfg_filename);
    return 0;
}
