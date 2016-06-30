#include <logger/logger.hpp>
#include <logger/log_writers.hpp>
#include <verify.hpp>
#include <segfault_handler.hpp>
#include <perfcounter.hpp>

#include <copy_file.hpp>

#include "vdj_config.hpp"
// #include "vdj_launch.hpp"

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

std::string running_time_format(const perf_counter &pc) {
    unsigned ms = (unsigned)pc.time_ms();
    unsigned secs = (ms / 1000) % 60;
    unsigned mins = (ms / 1000 / 60) % 60;
    unsigned hours = (ms / 1000 / 60 / 60);
    boost::format bf("%u hours %u minutes %u seconds");
    bf % hours % mins % secs;
    return bf.str();
}

void prepare_output_dir(const vdj_labeler::VDJLabelerConfig::IOParams::OutputParams & op) {
    path::make_dir(op.output_dir);
}

void copy_configs(std::string cfg_filename, std::string to) {
    path::make_dir(to);
    path::copy_files_by_ext(path::parent_path(cfg_filename), to, ".info", true);
    path::copy_files_by_ext(path::parent_path(cfg_filename), to, ".properties", true);
}

std::string get_config_fname(int argc, char **argv) {
    if(argc == 2 and (std::string(argv[1]) != "--help" and std::string(argv[1]) != "--version" and
                      std::string(argv[1]) != "--help-hidden"))
        return std::string(argv[1]);
    return "configs/vdj_labeler/config.info";
}


vdj_labeler::VDJLabelerConfig load_config(int argc, char **argv) {
    std::string cfg_filename = get_config_fname(argc, argv);
    vdj_labeler::VDJLabelerConfig vdj_config;
    vdj_config.load(cfg_filename);
    prepare_output_dir(vdj_config.io_params.output_params);
    std::string path_to_copy = path::append_path(vdj_config.io_params.output_params.output_dir, "configs");
    path::make_dir(path_to_copy);
    copy_configs(cfg_filename, path_to_copy);
    return vdj_config;
}


int main(int argc, char **argv) {
    segfault_handler sh;
    perf_counter pc;
    create_console_logger();
    vdj_labeler::VDJLabelerConfig config = load_config(argc, argv);
    // vdj_labeler::VDJLabelerLaunch(config).Launch();
    return 0;
}