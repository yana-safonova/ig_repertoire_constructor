#include <chrono>

#include "logger/log_writers.hpp"
#include "segfault_handler.hpp"
#include "copy_file.hpp"

#include "shm_kmer_matrix_estimator_config.hpp"
#include "shm_kmer_matrix_estimator_pipeline.hpp"

using namespace shm_kmer_matrix_estimator;

std::string get_config_fname(int argc, char **argv) {
    if(argc == 2)
        return std::string(argv[1]);
    return "configs/shm_kmer_matrix_estimator/config.info";
}

std::string load_config(int argc, char **argv) {
    std::string cfg_filename = get_config_fname(argc, argv);
    path::CheckFileExistenceFATAL(cfg_filename);
    shm_cfg::create_instance(cfg_filename);
    return cfg_filename;
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
        int error_code = shm_kmer_matrix_estimator::SHMkmerMatrixEstimatorPipeline(shm_cfg::get().io,
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
