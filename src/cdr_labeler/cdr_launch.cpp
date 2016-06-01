#include "cdr_launch.hpp"

#include <read_archive.hpp>
#include "germline_db_generator.hpp"
#include "germline_db_labeler.hpp"
#include "vj_parallel_processor.hpp"


namespace cdr_labeler {
    germline_utils::CustomGeneDatabase CDRLabelerLaunch::GetDatabaseByCDRLabeling(
            const germline_utils::CustomGeneDatabase &gene_db,
            DbCDRLabeling cdr_labeling) {
        germline_utils::CustomGeneDatabase labeled_db(gene_db.Segment());
        for(auto it = gene_db.cbegin(); it != gene_db.cend(); it++) {
            auto specific_gene_db = gene_db.GetDbByGeneType(*it);
            for(size_t i = 0; i < specific_gene_db.size(); i++) {
                auto gene_labeling = cdr_labeling.GetLabelingByGene(specific_gene_db[i]);
                if(!gene_labeling.Empty())
                    labeled_db.AddImmuneGene(specific_gene_db[i]);
            }
        }
        return labeled_db;
    }

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
        INFO("CDR labeling for germline segments");
        INFO("CDR labeling for V gene segments");
        auto v_labeling = GermlineDbLabeler(v_db, config_.cdrs_params).ComputeLabeling();
        INFO("CDR labeling for J gene segments");
        auto j_labeling = GermlineDbLabeler(j_db, config_.cdrs_params).ComputeLabeling();
        INFO("Creation of labeled V and J databases");
        auto labeled_v_db = GetDatabaseByCDRLabeling(v_db, v_labeling);
        auto labeled_j_db = GetDatabaseByCDRLabeling(j_db, j_labeling);
        INFO("Alignment against VJ germline segments");
        vj_finder::VJParallelProcessor processor(read_archive, config_.vj_finder_config.algorithm_params,
                                                 labeled_v_db, labeled_j_db,
                                                 config_.run_params.num_threads);
        vj_finder::VJAlignmentInfo alignment_info = processor.Process();
        INFO(alignment_info.NumVJHits() << " reads were aligned; " << alignment_info.NumFilteredReads() <<
             " reads were filtered out");
        INFO("CDR labeler ends");
    }
}