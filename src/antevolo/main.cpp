#include <logger/logger.hpp>
#include <logger/log_writers.hpp>
#include <verify.hpp>
#include <segfault_handler.hpp>
#include <perfcounter.hpp>

#include <copy_file.hpp>

#include "antevolo_config.hpp"
#include "antevolo_launch.hpp"
#include <vj_class_processors/edmonds_tarjan_DMST_calculator.hpp>

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

void prepare_output_dir(const antevolo::AntEvoloConfig &config) {
    path::make_dir(config.output_params.output_dir);
    path::make_dir(config.output_params.tree_dir);
    path::make_dir(config.output_params.vertex_dir);
    path::make_dir(config.output_params.cdr_graph_dir);
    path::make_dir(config.output_params.tree_shm_dir);
    if(config.algorithm_params.parallel_evolution_params.enable_parallel_shms_finder) {
        path::make_dir(config.output_params.parallel_shm_output.parallel_bulges_dir);
        path::make_dir(config.output_params.parallel_shm_output.parallel_shm_dir);
    }
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
    return "configs/antevolo/config.info";
}


antevolo::AntEvoloConfig load_config(int argc, char **argv) {
    std::string cfg_filename = get_config_fname(argc, argv);
    antevolo::AntEvoloConfig antevolo_config;
    antevolo_config.load(cfg_filename);
    prepare_output_dir(antevolo_config);
    std::string path_to_copy = path::append_path(antevolo_config.output_params.output_dir, "configs");
    path::make_dir(path_to_copy);
    copy_configs(cfg_filename, path_to_copy);
    return antevolo_config;
}


int main(int argc, char **argv) {
    omp_set_num_threads(1);

    segfault_handler sh;
    perf_counter pc;
    create_console_logger();

    antevolo::AntEvoloConfig config = load_config(argc, argv);
    antevolo::AntEvoloLaunch(config).Launch();

    return 0;
}
