#include "standard.hpp"
#include "logger/log_writers.hpp"

#include "segfault_handler.hpp"
#include "shm_config.hpp"
#include "copy_file.hpp"
#include "shm_kmer_model_estimator.hpp"

void make_dirs() {
    make_dir(shm_cfg::get().io.output.output_dir);
}

void copy_configs(string cfg_filename, string to) {
    if (!make_dir(to)) {
        WARN("Could not create files use in /tmp directory");
    }
    path::copy_files_by_ext(path::parent_path(cfg_filename), to, ".info", true);
}

void load_config(string cfg_filename) {
    path::CheckFileExistenceFATAL(cfg_filename);
    shm_cfg::create_instance(cfg_filename);
    string path_to_copy = path::append_path(shm_cfg::get().io.output.output_dir, "configs");
    copy_configs(cfg_filename, path_to_copy);
}

void create_console_logger(string cfg_filename) {
    using namespace logging;

    string log_props_file = shm_cfg::get().io.output.log_filename;

    if (!path::FileExists(log_props_file)){
        log_props_file = path::append_path(path::parent_path(cfg_filename),
                                           shm_cfg::get().io.output.log_filename);
    }

    logger *lg = create_logger(path::FileExists(log_props_file) ? log_props_file : "");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cout << "Invalid input parameters" << std::endl;
        std::cout << "shm_kmer_model config.info" << std::endl;
        return 1;
    }

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now() ;
    segfault_handler sh;

    try {
        string cfg_filename = argv[1];
        load_config(cfg_filename);
        create_console_logger(cfg_filename);
        make_dirs();
        int error_code =  shm_kmer_model_estimator::SHMkmerModelEstimator(shm_cfg::get().io,
                                                                          shm_cfg::get().achp,
                                                                          shm_cfg::get().acrp,
                                                                          shm_cfg::get().mfp).Run();
        if (error_code != 0) {
            INFO("SHM Kmer-Model Calculator finished abnormally");
            return error_code;
        }

    } catch (std::bad_alloc const& e) {
        std::cerr << "Not enough memory to run SHM Kmer-Model Calculator." << e.what() << std::endl;
        return EINTR;
    } catch (std::exception const& e) {
        std::cerr << "Exception caught " << e.what() << std::endl;
        return EINTR;
    } catch (...) {
        std::cerr << "Unknown exception caught " << std::endl;
        return EINTR;
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now() ;
    auto duration = end - start;
    auto secs = std::chrono::duration_cast<std::chrono::seconds>(duration);
    auto mins = std::chrono::duration_cast<std::chrono::minutes>(duration);
    auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
    INFO("Running time: " << hours.count() << " hours " << mins.count() <<
        " minutes " << secs.count() << " seconds");

    return 0;
}
