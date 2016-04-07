#include <logger/logger.hpp>

#include "vjf_launch.hpp"

#include <read_archive.hpp>
#include "germline_db_generator.hpp"


namespace vj_finder {
    void VJFinderLaunch::Run() {
        INFO("== VJ Finder starts == ");
        core::ReadArchive read_archive(config_.io_params.input_params.input_reads);
        GermlineDbGenerator db_generator(config_.io_params.input_params.germline_input,
                                         config_.algorithm_params.germline_params);
        INFO("Generation of DB for variable segments...");
        germline_utils::CustomGeneDatabase v_db = db_generator.GenerateVariableDb();
        INFO("Generation of DB for join segments...");
        germline_utils::CustomGeneDatabase j_db = db_generator.GenerateJoinDb();
        INFO("== VJ Finder ends == ");
    }
}