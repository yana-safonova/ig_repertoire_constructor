#include <chrono>

#include "logger/log_writers.hpp"
#include "segfault_handler.hpp"
#include "copy_file.hpp"

#include "shm_config.hpp"
#include "command_line_routines.hpp"
#include "shm_kmer_matrix_estimator.hpp"

using namespace shm_kmer_matrix_estimator;

void make_dirs() {
    path::make_dir(shm_cfg::get().io.output.output_dir);
}

void copy_configs(std::string cfg_filename, std::string to) {
    if (!path::make_dir(to)) {
        WARN("Could not create files use in /tmp directory");
    }
    path::copy_files_by_ext(path::parent_path(cfg_filename), to, ".info", true);
}

std::string get_config_fname(int argc, char **argv) {
    if(argc == 2 and (std::string(argv[1]) != "--help" and std::string(argv[1]) != "-h"))
        return std::string(argv[1]);
    return "configs/shm_kmer_model/configs.info";
}

std::string load_config(int argc, char **argv) {
    std::string cfg_filename = get_config_fname(argc, argv);
    path::CheckFileExistenceFATAL(cfg_filename);
    shm_cfg::create_instance(cfg_filename);
    std::string path_to_copy = path::append_path(shm_cfg::get().io.output.output_dir, "configs");
    path::make_dir(path_to_copy);
    copy_configs(cfg_filename, path_to_copy);
    parse_command_line_args(shm_cfg::get_writable(), argc, argv);
    return cfg_filename;
}

void load_config(std::string cfg_filename) {
    path::CheckFileExistenceFATAL(cfg_filename);
    shm_cfg::create_instance(cfg_filename);
    std::string path_to_copy = path::append_path(shm_cfg::get().io.output.output_dir, "configs");
    copy_configs(cfg_filename, path_to_copy);
}

void create_console_logger(std::string cfg_filename) {
    using namespace logging;

    std::string log_props_file = shm_cfg::get().io.output.log_filename;

    if (!path::FileExists(log_props_file)) {
        log_props_file = path::append_path(path::parent_path(cfg_filename),
                                           shm_cfg::get().io.output.log_filename);
    }

    logger *lg = create_logger(path::FileExists(log_props_file) ? log_props_file : "");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

int main(int argc, char *argv[]) {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    segfault_handler sh;

    try {
        std::string cfg_filename = load_config(argc, argv);
        create_console_logger(cfg_filename);
        int error_code = shm_kmer_matrix_estimator::SHMkmerMatrixEstimator(shm_cfg::get().io,
                                                                           shm_cfg::get().achp,
                                                                           shm_cfg::get().acrp,
                                                                           shm_cfg::get().mfp).Run();
        if (error_code != 0) {
            INFO("SHM k-mer Model Calculator finished abnormally. Error code: " << error_code);
            return error_code;
        }

    } catch (std::bad_alloc const &e) {
        std::cerr << "Not enough memory to run SHM Kmer-Model Calculator." << e.what() << std::endl;
        return EINTR;
    } catch (std::exception const &e) {
        std::cerr << "Exception caught " << e.what() << std::endl;
        return EINTR;
    } catch (...) {
        std::cerr << "Unknown exception caught " << std::endl;
        return EINTR;
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    auto duration = end - start;
    auto secs = std::chrono::duration_cast<std::chrono::seconds>(duration);
    auto mins = std::chrono::duration_cast<std::chrono::minutes>(duration);
    auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
    INFO("Running time: " << hours.count() << " hours " << mins.count() <<
        " minutes " << secs.count() << " seconds");

    return 0;
}
