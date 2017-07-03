#include <logger/logger.hpp>

#include "vjf_launch.hpp"

#include <read_archive.hpp>
#include "germline_utils/germline_db_generator.hpp"
#include "vj_parallel_processor.hpp"

using namespace germline_utils;

namespace vj_finder {
    void CreateAlignmentOutput(std::ofstream& fhandler, const core::Read& read, const VJHits& vj_hits) {
        fhandler << read.name << "\t" << vj_hits.GetVHitByIndex(0).Start() << "\t" <<
                vj_hits.GetVHitByIndex(0).End() << "\t" << vj_hits.GetVHitByIndex(0).Score() << "\t" <<
                vj_hits.GetVHitByIndex(0).ImmuneGene().name() << "\t" << vj_hits.GetJHitByIndex(0).Start() << "\t" <<
                vj_hits.GetJHitByIndex(0).End() << "\t" << vj_hits.GetJHitByIndex(0).Score() << "\t" <<
                vj_hits.GetJHitByIndex(0).ImmuneGene().name() << std::endl;
    }

    void VJFinderLaunch::Run() {
        INFO("== VJ Finder starts == ");
        core::ReadArchive read_archive(config_.io_params.input_params.input_reads);
        if(config_.io_params.output_params.output_details.fix_spaces)
            read_archive.FixSpacesInHeaders();
        GermlineDbGenerator db_generator(config_.io_params.input_params.germline_input,
                                         config_.algorithm_params.germline_params);
        INFO("Generation of DB for variable segments...");
        germline_utils::CustomGeneDatabase v_db = db_generator.GenerateVariableDb();
        INFO("Generation of DB for join segments...");
        germline_utils::CustomGeneDatabase j_db = db_generator.GenerateJoinDb();
        VJParallelProcessor processor(read_archive, config_.algorithm_params, v_db, j_db,
                                      config_.run_params.num_threads);
        INFO("Alignment against VJ germline segments starts");
        VJAlignmentInfo alignment_info = processor.Process();
        INFO(alignment_info.NumVJHits() << " reads were aligned; " << alignment_info.NumFilteredReads() <<
                     " reads were filtered out");
        VJAlignmentOutput alignment_info_output(config_.io_params.output_params, alignment_info);
        alignment_info_output.OutputAlignmentInfo();
        alignment_info_output.OutputCleanedReads();
        alignment_info_output.OutputVAlignments();
        alignment_info_output.OutputFilteredReads();
        alignment_info_output.OutputFilteringInfo();
        INFO("== VJ Finder ends == ");
    }
}