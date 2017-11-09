#include <logger/logger.hpp>
#include <logger/log_writers.hpp>
#include <segfault_handler.hpp>

#include <copy_file.hpp>

#include "ig_simulator_config.hpp"
#include "ig_simulator_launch.hpp"
#include "../config/config_utils.hpp"

namespace {
    void create_console_logger(std::string cfg_filename) {
        using namespace logging;
        std::string log_props_file = ig_simulator::igs_cfg::get().io_params.output_params.log_filename;
        if (!path::FileExists(log_props_file)){
            log_props_file = path::append_path(path::parent_path(cfg_filename), log_props_file);
        }
        logger *lg = create_logger(path::FileExists(log_props_file) ? log_props_file : "");
        lg->add_writer(std::make_shared<console_writer>());
        attach_logger(lg);
    }

    std::string running_time_format(const perf_counter &pc) {
        unsigned ms = (unsigned)pc.time_ms();
        unsigned secs = (ms / 1000) % 60;
        unsigned mins = (ms / 1000 / 60) % 60;
        unsigned hours = (ms / 1000 / 60 / 60);
        boost::format bf("%u hours %u minutes %u seconds");
        bf % hours % mins % secs;
        return bf.str();
    }

    class IgsConfigLoader : public config_utils::ConfigLoader<ig_simulator::IgSimulatorConfig> {
    protected:
        void FillConfigFromCommandline(ig_simulator::IgSimulatorConfig&, int, const char* const*) const override {}

        std::string GetDefaultCfgFilename() const override {
            return "configs/ig_simulator/config.info";
        }

        std::string GetOutputDirPath(const ig_simulator::IgSimulatorConfig& config) const override {
            return config.io_params.output_params.output_dir;
        }
    };
}

int main(int argc, char **argv) {
    omp_set_num_threads(1);

    segfault_handler sh;
    perf_counter pc;
    std::string cfg_filename = IgsConfigLoader().LoadConfig(argc, argv);
    create_console_logger(cfg_filename);
    // variable extracted to avoid a possible bug in gcc 4.8.4
    const auto& cfg = ig_simulator::igs_cfg::get();
    ig_simulator::IgSimulatorLaunch(cfg).Run();
    return 0;
}
