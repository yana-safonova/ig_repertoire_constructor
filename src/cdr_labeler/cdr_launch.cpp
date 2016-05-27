#include "cdr_launch.hpp"

#include <read_archive.hpp>
#include "germline_db_generator.hpp"
#include "germline_db_labeler.hpp"
#include "vj_parallel_processor.hpp"


namespace cdr_labeler {
    void CDRLabelerLaunch::Launch() {
        INFO("CDR labeler starts");
        core::ReadArchive read_archive(config_.input_params.input_reads);
        if(config_.vj_finder_config.io_params.output_params.output_details.fix_spaces)
            read_archive.FixSpacesInHeaders();
        vj_finder::GermlineDbGenerator db_generator(config_.vj_finder_config.io_params.input_params.germline_input,
                                         config_.vj_finder_config.algorithm_params.germline_params);
        INFO("Generation of DB for variable segments...");
        germline_utils::CustomGeneDatabase v_db = db_generator.GenerateVariableDb();
        INFO("Generation of DB for join segments...");
        germline_utils::CustomGeneDatabase j_db = db_generator.GenerateJoinDb();
        INFO("CDR labeling for germline segments starts");
        auto v_labeling = GermlineDbLabeler(v_db, config_.cdrs_params).ComputeLabeling();
        INFO("Alignment against VJ germline segments starts");
        vj_finder::VJParallelProcessor processor(read_archive, config_.vj_finder_config.algorithm_params, v_db, j_db,
                                                 config_.run_params.num_threads);
        vj_finder::VJAlignmentInfo alignment_info = processor.Process();
        INFO(alignment_info.NumVJHits() << " reads were aligned; " << alignment_info.NumFilteredReads() <<
             " reads were filtered out");
        INFO("CDR labeler ends");
    }
}