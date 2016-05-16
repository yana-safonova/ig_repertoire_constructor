#include <logger/logger.hpp>

#include "vjf_launch.hpp"

#include <read_archive.hpp>
#include "germline_db_generator.hpp"
#include "vj_query_aligner.hpp"
#include "vj_hits_filter.hpp"


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
        VJQueryAligner vj_query_aligner(config_.algorithm_params, v_db, j_db);
        VersatileVjFilter vj_filter(config_.algorithm_params.filtering_params);
        std::string alignment_info = config_.io_params.output_params.output_files.add_info_filename;
        std::ofstream alignment_info_fhandler(alignment_info);
        for(auto it = read_archive.cbegin(); it != read_archive.cend(); it++) {
            TRACE("Processing read: " << it->name << ", id: " << it->id << ", length: " << it->length());
            auto vj_hits = vj_query_aligner.Align(*it);
            if(!vj_filter.Filter(vj_hits)) {
                TRACE("Read is good");
                CreateAlignmentOutput(alignment_info_fhandler, *it, vj_hits);
            }
            else
                TRACE("Read was filtered out");
        }
        INFO("Alignment information was written to " << alignment_info);
        alignment_info_fhandler.close();
        INFO("== VJ Finder ends == ");
    }
}