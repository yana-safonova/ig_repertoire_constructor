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

std::string running_time_format(const perf_counter &pc) {
    unsigned ms = (unsigned)pc.time_ms();
    unsigned secs = (ms / 1000) % 60;
    unsigned mins = (ms / 1000 / 60) % 60;
    unsigned hours = (ms / 1000 / 60 / 60);
    boost::format bf("%u hours %u minutes %u seconds");
    bf % hours % mins % secs;
    return bf.str();
}

void prepare_output_dir(const vj_finder::VJFinderConfig::IOParams::OutputParams::OutputFiles & of) {
    path::make_dir(of.output_dir);
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
    return "configs/vj_finder/config.info";
}

std::string load_config(int argc, char **argv) {
    std::string cfg_filename = get_config_fname(argc, argv);
    if (!path::FileExists(cfg_filename)) {
        std::cout << "File " << cfg_filename << " doesn't exist or can't be read!" << std::endl;
        exit(-1);
    }
    vj_finder::vjf_cfg::create_instance(cfg_filename);
    parse_command_line_args(vj_finder::vjf_cfg::get_writable(), argc, argv);
    prepare_output_dir(vj_finder::vjf_cfg::get().io_params.output_params.output_files);
    std::string path_to_copy = path::append_path(vj_finder::vjf_cfg::get().io_params.output_params.output_files.output_dir, "configs");
    path::make_dir(path_to_copy);
    copy_configs(cfg_filename, path_to_copy);
    return cfg_filename;
}

int main(int argc, char **argv) {
    segfault_handler sh;
    perf_counter pc;
    std::string cfg_filename = load_config(argc, argv);
    create_console_logger(cfg_filename);
    // variable extracted to avoid a possible bug in gcc 4.8.4
    const auto& cfg = vj_finder::vjf_cfg::get();
    vj_finder::VJFinderLaunch(cfg).Run();
    return 0;
}