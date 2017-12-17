#include "logger/log_writers.hpp"
#include "omp.h"

#include "segfault_handler.hpp"
#include "stacktrace.hpp"
#include "memory_limit.hpp"
#include "copy_file.hpp"
#include "perfcounter.hpp"
#include "runtime_k.hpp"
#include "segfault_handler.hpp"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "path_helper.hpp"
#include "dsf_config.hpp"
#include "launch.hpp"
#include "../config/config_utils.hpp"

namespace {
    void make_dirs() {
        using path::make_dir;
        make_dir(dsf_cfg::get().io.output_base.output_dir);
        if(dsf_cfg::get().rp.threads_count > 1) {
            make_dir(dsf_cfg::get().io.output_mthreading.connected_components_dir);
            make_dir(dsf_cfg::get().io.output_mthreading.decompositions_dir);
        }
    }

    void create_console_logger(std::string cfg_filename) {
        using namespace logging;

        std::string log_props_file = dsf_cfg::get().io.output_base.log_filename;

        if (!path::FileExists(log_props_file)){
            log_props_file = path::append_path(path::parent_path(cfg_filename), dsf_cfg::get().io.output_base.log_filename);
        }

        logger *lg = create_logger(path::FileExists(log_props_file) ? log_props_file : "");
        lg->add_writer(std::make_shared<console_writer>());
        attach_logger(lg);
    }

    class DsfConfigLoader : public config_utils::ConfigLoader<dsf_config> {
    protected:
        void FillConfigFromCommandline(dsf_config&, int, const char* const*) const override {}

        void PrepareOutputDir(const std::string&) const override {}

        std::string GetDefaultCfgFilename() const override {
            std::cout << "Usage: dense_sgraph_finder <config file>\n";
            exit(-1);
        }

        std::string GetOutputDirPath(const dsf_config& config) const override {
            return config.io.output_base.output_dir;
        }

        void CopyConfigs(const std::string& cfg_filename, const dsf_config& config) const override {
            //using namespace debruijn_graph;

            std::string to_dir = path::append_path(GetOutputDirPath(config), "configs");
            if (!path::make_dir(to_dir)) {
                WARN("Could not create files use in /tmp directory");
            }
            path::copy_files_by_ext(path::parent_path(cfg_filename), to_dir, ".info", true);
        }
    };
}

int main(int argc, char* argv[]) {
    omp_set_num_threads(1);

    if(argc != 2) {
        std::cout << "dense_sgraph_finder config.info" << std::endl;
        return 1;
    }

    perf_counter pc;
    segfault_handler sh;

    try {
        std::string cfg_filename = DsfConfigLoader().LoadConfig(argc, argv);
        create_console_logger(cfg_filename);
        make_dirs();
        int error_code = dense_subgraph_finder::DenseSubgraphFinder(dsf_cfg::get().rp,
                                                                    dsf_cfg::get().dsf_params,
                                                                    dsf_cfg::get().io,
                                                                    dsf_cfg::get().metis_io).Run();
        if (error_code != 0) {
            INFO("Dense subgraph finder finished abnormally");
            return error_code;
        }
    } catch (std::bad_alloc const& e) {
        std::cerr << "Not enough memory to run IgRepertoireConstructor. " << e.what() << std::endl;
        return EINTR;
    } catch (std::exception const& e) {
        std::cerr << "Exception caught " << e.what() << std::endl;
        return EINTR;
    } catch (...) {
        std::cerr << "Unknown exception caught " << std::endl;
        return EINTR;
    }

    unsigned ms = (unsigned)pc.time_ms();
    unsigned secs = (ms / 1000) % 60;
    unsigned mins = (ms / 1000 / 60) % 60;
    unsigned hours = (ms / 1000 / 60 / 60);
    INFO("Running time: " << hours << " hours " << mins << " minutes " << secs << " seconds");

    return 0;
}
