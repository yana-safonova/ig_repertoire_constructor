#include <logger/logger.hpp>
#include <logger/log_writers.hpp>
#include <segfault_handler.hpp>

#include <copy_file.hpp>

#include "ig_simulator_config.hpp"
#include "ig_simulator_launch.hpp"

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

void prepare_output_dir(const ig_simulator::IgSimulatorConfig::IOParams::OutputParams & of) {
    path::make_dir(of.output_dir);
}

void copy_configs(std::string cfg_filename, std::string to) {
    path::make_dir(to);
    path::copy_files_by_ext(path::parent_path(cfg_filename), to, ".info", true);
    path::copy_files_by_ext(path::parent_path(cfg_filename), to, ".properties", true);
}

std::string get_config_fname(int argc, char **argv) {
    if(argc == 2)
        return std::string(argv[1]);
    return "configs/ig_simulator/config.info";
}

std::string load_config(int argc, char **argv) {
    std::string cfg_filename = get_config_fname(argc, argv);
    if (!path::FileExists(cfg_filename)) {
        std::cout << "File " << cfg_filename << " doesn't exist or can't be read!" << std::endl;
        exit(-1);
    }
    ig_simulator::igs_cfg::create_instance(cfg_filename);
    prepare_output_dir(ig_simulator::igs_cfg::get().io_params.output_params);
    std::string path_to_copy =
        path::append_path(ig_simulator::igs_cfg::get().io_params.output_params.output_dir, "configs");
    path::make_dir(path_to_copy);
    copy_configs(cfg_filename, path_to_copy);
    return cfg_filename;
}

int main(int argc, char **argv) {
    omp_set_num_threads(1);

    segfault_handler sh;
    perf_counter pc;
    std::string cfg_filename = load_config(argc, argv);
    create_console_logger(cfg_filename);
    // variable extracted to avoid a possible bug in gcc 4.8.4
    const auto& cfg = ig_simulator::igs_cfg::get();
    ig_simulator::IgSimulatorLaunch(cfg).Run();
    return 0;
}
