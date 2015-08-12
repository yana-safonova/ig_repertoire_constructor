#include "standard.hpp"
#include "logger/log_writers.hpp"

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

#include "dsf_config.hpp"
#include "launch.hpp"

void make_dirs() {
    make_dir(dsf_cfg::get().io.output_base.output_dir);
    if(dsf_cfg::get().rp.threads_count > 1) {
        make_dir(dsf_cfg::get().io.output_mthreading.connected_components_dir);
        make_dir(dsf_cfg::get().io.output_mthreading.decompositions_dir);
    }
}

void copy_configs(string cfg_filename, string to) {
    //using namespace debruijn_graph;

    if (!make_dir(to)) {
        WARN("Could not create files use in /tmp directory");
    }
    path::copy_files_by_ext(path::parent_path(cfg_filename), to, ".info", true);
}

void load_config(string cfg_filename) {
    path::CheckFileExistenceFATAL(cfg_filename);
    dsf_cfg::create_instance(cfg_filename);
    string path_to_copy = path::append_path(dsf_cfg::get().io.output_base.output_dir, "configs");
    copy_configs(cfg_filename, path_to_copy);
}

void create_console_logger(string cfg_filename) {
    using namespace logging;

    string log_props_file = dsf_cfg::get().io.output_base.log_filename;

    if (!path::FileExists(log_props_file)){
        log_props_file = path::append_path(path::parent_path(cfg_filename), dsf_cfg::get().io.output_base.log_filename);
    }

    logger *lg = create_logger(path::FileExists(log_props_file) ? log_props_file : "");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

int main(int argc, char* argv[]) {
    if(argc != 2) {
        std::cout << "dense_sgraph_finder config.info" << std::endl;
        return 1;
    }

    string cfg_filename = argv[1];
    load_config(cfg_filename);
    create_console_logger(cfg_filename);
    make_dirs();

    int error_code = dense_subgraph_finder::DenseSubgraphFinder().Run(dsf_cfg::get().dsf_params,
                                                                      dsf_cfg::get().io,
                                                                      dsf_cfg::get().metis_io);
    if(error_code != 0) {
        INFO("Dense subgraph finder finished abnormally");
        return error_code;
    }
    return 0;
}