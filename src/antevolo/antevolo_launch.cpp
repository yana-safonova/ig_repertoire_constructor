#include "antevolo_launch.hpp"

#include <read_archive.hpp>
#include <germline_db_generator.hpp>
#include <cdr_annotator.hpp>

#include "naive_antevolo_processing.hpp"

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
        auto annotated_clone_set = cdr_labeler::CDRAnnotator(config_.cdr_labeler_config,
                                                             read_archive, v_db, j_db).AnnotateClones();
        INFO("Naive tree construction starts");
        NaiveAntEvoloProcessing(config_, annotated_clone_set).ConstructClonalTrees();
        INFO("AntEvolo ends");
    }
}