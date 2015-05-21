//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "standard.hpp"
#include "logger/log_writers.hpp"

#include "segfault_handler.hpp"
#include "stacktrace.hpp"
#include "memory_limit.hpp"
#include "copy_file.hpp"
#include "perfcounter.hpp"
#include "runtime_k.hpp"
#include "io/ireadstream.hpp"
#include "io/library.hpp"
#include "io/single_read.hpp"
#include "io/io_helper.hpp"

//#include "config_struct.hpp"

#include "stage.hpp"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "segfault_handler.hpp"

#include "ig_config.hpp"
#include "launch.hpp"

void make_dirs(){
    make_dir(ig_cfg::get().io.output_dir);
    make_dir(ig_cfg::get().io.temp_files);
    if (ig_cfg::get().rp.developer_mode)
      make_dir(ig_cfg::get().io.output_saves);
    make_dir(ig_cfg::get().hgc_params.hgc_io_params.hg_output_dir);
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
    ig_cfg::create_instance(cfg_filename);
    string path_to_copy = path::append_path(ig_cfg::get().io.output_dir, "configs");
    copy_configs(cfg_filename, path_to_copy);
}

void create_console_logger(string cfg_filename) {
  using namespace logging;

  string log_props_file = ig_cfg::get().io.log_filename;

  if (!path::FileExists(log_props_file)){
    log_props_file = path::append_path(path::parent_path(cfg_filename), ig_cfg::get().io.log_filename);
  }

  logger *lg = create_logger(path::FileExists(log_props_file) ? log_props_file : "");
  lg->add_writer(std::make_shared<console_writer>());
  attach_logger(lg);
}

int main(int /*argc*/, char** argv) {
    perf_counter pc;

    srand(42);
    srandom(42);

    segfault_handler sh;

    try {
        //using namespace debruijn_graph;
        string cfg_filename = argv[1];
        load_config          (cfg_filename);
        make_dirs();
        copy_configs(cfg_filename, path::append_path(ig_cfg::get().io.output_dir, "configs"));
        create_console_logger(cfg_filename);

        const size_t GB = 1 << 30;
        limit_memory(ig_cfg::get().rp.max_memory * GB);

        ig_repertoire_constructor::run_ig_repertoire_constructor();

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

    // OK
    return 0;
}
