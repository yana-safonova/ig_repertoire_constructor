#include <logger/logger.hpp>

#include "vjf_launch.hpp"

#include <read_archive.hpp>
#include "germline_db_generator.hpp"
#include "vj_query_aligner.hpp"
#include "vj_hits_filter.hpp"


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
        VJQueryAligner vj_query_aligner(config_.algorithm_params, v_db, j_db);
        VersatileVjFilter vj_filter(config_.algorithm_params.filtering_params);
        for(auto it = read_archive.cbegin(); it != read_archive.cend(); it++) {
            std::cout << "Read: " << it->name << ", id: " << it->id << ", length: " << it->length() << std::endl;
            auto vj_hits = vj_query_aligner.Align(*it);
            if(!vj_filter.Filter(vj_hits))
                std::cout << "Read is good" << std::endl;
            else
                std::cout << "Read to be filtered out" << std::endl;
            std::cout << "------------------------" << std::endl;
        }
        INFO("== VJ Finder ends == ");
    }
}