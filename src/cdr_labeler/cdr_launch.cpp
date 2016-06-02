#include "cdr_launch.hpp"

#include <read_archive.hpp>
#include "germline_db_generator.hpp"
#include "germline_db_labeler.hpp"
#include "vj_parallel_processor.hpp"
#include "immune_gene_alignment_converter.hpp"

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
        INFO("CDR labeling for V gene segments");
        auto v_labeling = GermlineDbLabeler(v_db, config_.cdrs_params).ComputeLabeling();
        INFO("CDR labeling for J gene segments");
        auto j_labeling = GermlineDbLabeler(j_db, config_.cdrs_params).ComputeLabeling();
        INFO("Creation of labeled V and J databases");
        auto labeled_v_db = v_labeling.CreateFilteredDb();
        INFO("Labeled DB of V segments consists of " << labeled_v_db.size() << " records");
        auto labeled_j_db = j_labeling.CreateFilteredDb();
        INFO("Labeled DB of J segments consists of " << labeled_j_db.size() << " records");
        INFO("Alignment against VJ germline segments");
        vj_finder::VJParallelProcessor processor(read_archive, config_.vj_finder_config.algorithm_params,
                                                 labeled_v_db, labeled_j_db,
                                                 config_.run_params.num_threads);
        vj_finder::VJAlignmentInfo alignment_info = processor.Process();
        INFO(alignment_info.NumVJHits() << " reads were aligned; " << alignment_info.NumFilteredReads() <<
                     " reads were filtered out");

        for(size_t i = 0; i < alignment_info.NumVJHits(); i++) {
            auto vj_hit = alignment_info.GetVJHitsByIndex(i);
            auto v_hit = vj_hit.GetVHitByIndex(0);
            vj_finder::ImmuneGeneAlignmentConverter converter;
            auto v_alignment = converter.ConvertToAlignment(v_hit.ImmuneGene(), v_hit.Read(), v_hit.BlockAlignment());
            std::cout << v_hit.ImmuneGene() << std::endl;
            std::cout << v_hit.Read() << std::endl;
            std::cout << v_alignment.Alignment() << std::endl;
            auto v_cdr_labeling = v_labeling.GetLabelingByGene(v_hit.ImmuneGene());
            CDRRange read_cdr1(v_alignment.GetQueryPositionBySubjectPosition(v_cdr_labeling.cdr1.start_pos),
                               v_alignment.GetQueryPositionBySubjectPosition(v_cdr_labeling.cdr1.end_pos));
            std::cout << seqan::infixWithLength(v_hit.Read().seq, read_cdr1.start_pos, read_cdr1.length()) << std::endl;

            auto j_hit = vj_hit.GetJHitByIndex(0);
            auto j_alignment = converter.ConvertToAlignment(j_hit.ImmuneGene(), j_hit.Read(), j_hit.BlockAlignment());
            std::cout << j_hit.ImmuneGene() << std::endl;
            std::cout << j_hit.Read() << std::endl;
            std::cout << j_alignment.Alignment() << std::endl;

            auto j_cdr_labeling = j_labeling.GetLabelingByGene(j_hit.ImmuneGene());
            CDRRange read_cdr3(v_alignment.GetQueryPositionBySubjectPosition(v_cdr_labeling.cdr3.start_pos),
                               j_alignment.GetQueryPositionBySubjectPosition(j_cdr_labeling.cdr3.end_pos));
            std::cout << seqan::infixWithLength(v_hit.Read().seq, read_cdr3.start_pos, read_cdr3.length()) << std::endl;
        }


        INFO("CDR labeler ends");
    }
}