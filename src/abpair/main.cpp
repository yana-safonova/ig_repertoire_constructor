#include "standard.hpp"
#include "logger/log_writers.hpp"

#include "segfault_handler.hpp"
#include "stacktrace.hpp"
#include "memory_limit.hpp"
#include "copy_file.hpp"
#include "perfcounter.hpp"
#include "runtime_k.hpp"
#include "segfault_handler.hpp"

#include "abpair_config.hpp"
#include "abpair_launch.hpp"

void make_dirs() {
    make_dir(abp_cfg::get().io.output.output_dir);
    if(abp_cfg::get().io.output.output_statistics) {
        make_dir(abp_cfg::get().io.output.hc_ambiguous_dir);
    }
    if(abp_cfg::get().io.output.output_barcodes) {
        make_dir(abp_cfg::get().io.output.barcode_dir);
    }
    if(abp_cfg::get().io.output.output_demultiplexed_raw_data)
        make_dir(abp_cfg::get().io.output.demultiplexed_raw_dir);
}

void copy_configs(std::string cfg_filename, std::string to) {
    if (!make_dir(to)) {
        WARN("Could not create files use in /tmp directory");
    }
    path::copy_files_by_ext(path::parent_path(cfg_filename), to, ".info", true);
}

void load_config(std::string cfg_filename) {
    path::CheckFileExistenceFATAL(cfg_filename);
    abp_cfg::create_instance(cfg_filename);
    std::string path_to_copy = path::append_path(abp_cfg::get().io.output.output_dir, "configs");
    copy_configs(cfg_filename, path_to_copy);
}

void create_console_logger(std::string cfg_filename) {
    using namespace logging;
    std::string log_props_file = abp_cfg::get().io.output.log_filename;
    if (!path::FileExists(log_props_file)){
        log_props_file = path::append_path(path::parent_path(cfg_filename), abp_cfg::get().io.output.log_filename);
    }
    logger *lg = create_logger(path::FileExists(log_props_file) ? log_props_file : "");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

int main(int argc, char** argv) {
    if(argc != 2) {
        std::cout << "Invalid input parameters" << std::endl;
        std::cout << "abpair config.info" << std::endl;
        return 1;
    }

    perf_counter pc;
    segfault_handler sh;

    std::string cfg_filename = argv[1];
    load_config(cfg_filename);
    create_console_logger(cfg_filename);
    make_dirs();

    AbPairLauncher().Run(abp_cfg::get().io);

    unsigned ms = (unsigned)pc.time_ms();
    unsigned secs = (ms / 1000) % 60;
    unsigned mins = (ms / 1000 / 60) % 60;
    unsigned hours = (ms / 1000 / 60 / 60);
    INFO("Running time: " << hours << " hours " << mins << " minutes " << secs << " seconds");

    return 0;
}