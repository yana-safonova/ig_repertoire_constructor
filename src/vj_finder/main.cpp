#include <logger/logger.hpp>
#include <logger/log_writers.hpp>
#include <verify.hpp>
#include <segfault_handler.hpp>
#include <perfcounter.hpp>

#include <read_archive.hpp>
#include <copy_file.hpp>

#include "vj_finder_config.hpp"
#include "command_line_routines.hpp"
#include "vjf_launch.hpp"
#include "../config/config_utils.hpp"

namespace {
    void create_console_logger(std::string cfg_filename) {
        using namespace logging;
        std::string log_props_file = vj_finder::vjf_cfg::get().io_params.output_params.output_files.log_filename;
        if (!path::FileExists(log_props_file)){
            log_props_file = path::append_path(path::parent_path(cfg_filename), log_props_file);
        }
        logger *lg = create_logger(path::FileExists(log_props_file) ? log_props_file : "");
        lg->add_writer(std::make_shared<console_writer>());
        attach_logger(lg);
    }

    class VjfConfigLoader : public config_utils::ConfigLoader<vj_finder::VJFinderConfig> {
    protected:
        std::string GetDefaultCfgFilename() const override {
            return "configs/vj_finder/config.info";
        }

        void FillConfigFromCommandline(vj_finder::VJFinderConfig& config, int argc, const char* const* argv) const override {
            parse_command_line_args(config, argc, argv);
        }

        std::string GetOutputDirPath(const vj_finder::VJFinderConfig& config) const override {
            return config.io_params.output_params.output_files.output_dir;
        }
    };
}


int main(int argc, const char* const* argv) {
    omp_set_num_threads(1);

    segfault_handler sh;
    perf_counter pc;
    std::string cfg_filename = VjfConfigLoader().LoadConfig(argc, argv);
    create_console_logger(cfg_filename);
    // variable extracted to avoid a possible bug in gcc 4.8.4
    const auto& cfg = vj_finder::vjf_cfg::get();
    vj_finder::VJFinderLaunch(cfg).Run();
    return 0;
}
