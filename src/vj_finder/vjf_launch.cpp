#include <logger/logger.hpp>

#include "vjf_launch.hpp"

#include <read_archive.hpp>
#include "germline_db_generator.hpp"
#include "vj_alignment_info.hpp"
#include "vj_query_processing.hpp"

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
        VJAlignmentInfo alignment_info;
        for(auto it = read_archive.cbegin(); it != read_archive.cend(); it++) {
            TRACE("Processing read: " << it->name << ", id: " << it->id << ", length: " << it->length());
            VJQueryProcessor vj_query_processor(config_.algorithm_params, v_db, j_db);
            auto vj_hits = vj_query_processor.Process(*it);
            if(vj_hits)
                alignment_info.UpdateFilteredReads(*it);
            else
                alignment_info.UpdateHits(*vj_hits);
        }
        VJAlignmentOutput alignment_info_output(config_.io_params.output_params, alignment_info);
        alignment_info_output.OutputAlignmentInfo();
        alignment_info_output.OutputCleanedReads();
        alignment_info_output.OutputFilteredReads();
        INFO("== VJ Finder ends == ");
    }
}