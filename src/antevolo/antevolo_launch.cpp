#include "antevolo_launch.hpp"

#include <read_archive.hpp>
#include <germline_db_generator.hpp>
#include <germline_db_labeler.hpp>
#include <vj_parallel_processor.hpp>
#include <read_labeler.hpp>
#include <cdr_output.hpp>

#include "antevolo_processor.hpp"
#include "evolutionary_stats_calculator.hpp"
#include "antevolo_output_writer.hpp"

namespace antevolo {
    void AntEvoloLaunch::Launch() {
        INFO("AntEvolo starts");
        core::ReadArchive read_archive(config_.input_params.input_reads);
        if(config_.cdr_labeler_config.vj_finder_config.io_params.output_params.output_details.fix_spaces)
            read_archive.FixSpacesInHeaders();
        vj_finder::GermlineDbGenerator db_generator(config_.cdr_labeler_config.vj_finder_config.io_params.input_params.germline_input,
                                                    config_.cdr_labeler_config.vj_finder_config.algorithm_params.germline_params);
        INFO("Generation of DB for variable segments...");
        germline_utils::CustomGeneDatabase v_db = db_generator.GenerateVariableDb();
        INFO("Generation of DB for join segments...");
        germline_utils::CustomGeneDatabase j_db = db_generator.GenerateJoinDb();
        // todo: refactor code duplication
        INFO("CDR labeling for V gene segments");
        auto v_labeling = cdr_labeler::GermlineDbLabeler(v_db,
                                                         config_.cdr_labeler_config.cdrs_params).ComputeLabeling();
        INFO("CDR labeling for J gene segments");
        auto j_labeling = cdr_labeler::GermlineDbLabeler(j_db,
                                                         config_.cdr_labeler_config.cdrs_params).ComputeLabeling();
        INFO("Creation of labeled V and J databases");
        auto labeled_v_db = v_labeling.CreateFilteredDb();
        INFO("Labeled DB of V segments consists of " << labeled_v_db.size() << " records");
        auto labeled_j_db = j_labeling.CreateFilteredDb();
        INFO("Labeled DB of J segments consists of " << labeled_j_db.size() << " records");
        INFO("Alignment against VJ germline segments");
        vj_finder::VJParallelProcessor processor(read_archive,
                                                 config_.cdr_labeler_config.vj_finder_config.algorithm_params,
                                                 labeled_v_db, labeled_j_db,
                                                 config_.cdr_labeler_config.run_params.num_threads);
        vj_finder::VJAlignmentInfo alignment_info = processor.Process();
        INFO(alignment_info.NumVJHits() << " reads were aligned; " << alignment_info.NumFilteredReads() <<
             " reads were filtered out");

        cdr_labeler::ReadCDRLabeler read_labeler(config_.cdr_labeler_config.shm_params, v_labeling, j_labeling);
        auto annotated_clone_set = read_labeler.CreateAnnotatedCloneSet(alignment_info);
        cdr_labeler::CDRLabelingWriter writer(config_.cdr_labeler_config.output_params,
                                              annotated_clone_set);
        writer.OutputCDRDetails();
        writer.OutputSHMs();
        INFO("Tree construction starts");
        auto tree_storage = AntEvoloProcessor(config_, annotated_clone_set).ConstructClonalTrees();
        INFO(tree_storage.size() << " evolutionary trees were created");
        INFO("Computation of evolutionary statistics");
        // todo: implement splitter into connected components
        EvolutionaryStatsCalculator stats_calculator(annotated_clone_set);
        auto annotated_storage = stats_calculator.ComputeStatsForStorage(tree_storage);
        INFO(annotated_storage.size() << " trees were annotated");
        AntEvoloOutputWriter output_writer(config_.output_params, annotated_storage);
        output_writer.OutputTreeStats();
        output_writer.OutputSHMForTrees();
        INFO("AntEvolo ends");
    }
}
