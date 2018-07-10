#include "cdr_launch.hpp"

#include <read_archive.hpp>
#include "germline_utils/germline_db_generator.hpp"
#include "germline_db_labeler.hpp"
#include "vj_parallel_processor.hpp"
#include "read_labeler.hpp"
#include "cdr_output.hpp"
#include "diversity_analyser.hpp"
#include "germline_utils/germline_config.hpp"
//#include "cdr_annotator.hpp"

namespace cdr_labeler {
    void CDRLabelerLaunch::Launch() {
        using namespace annotation_utils;
        CheckInputParams();
        INFO("Diversity Analyzer starts");
        core::ReadArchive read_archive(config_.input_params.input_reads);
        if(config_.vj_finder_config.io_params.output_params.output_details.fix_spaces)
            read_archive.FixSpacesInHeaders();
        germline_utils::GermlineDbGenerator db_generator(config_.vj_finder_config.io_params.input_params.germline_input,
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
        ReadCDRLabeler read_labeler(config_.shm_params, v_labeling, j_labeling);
        auto annotated_clone_set = read_labeler.CreateAnnotatedCloneSet(alignment_info);
        INFO("CDR sequences and SHMs were computed");
        CDRLabelingWriter writer(config_.output_params, annotated_clone_set);
        writer.OutputCleanedReads();
        writer.OutputCDRDetails();
        writer.OutputCDR1Fasta();
        writer.OutputCDR2Fasta();
        writer.OutputCDR3Fasta();
        writer.OutputCompressedCDR3Fasta();
        writer.OutputVGeneAlignment();
        writer.OutputSHMs();
        INFO("Diversity analysis of CDRs");
        DiversityAnalyser cdr_analyser(annotated_clone_set, config_.input_params,
                                       config_.output_params,
                                       config_.output_params.cdr3_compressed_fasta);
        INFO("Shannon index. CDR1: " << cdr_analyser.ShannonIndex(StructuralRegion::CDR1) <<
                ", CDR2: " << cdr_analyser.ShannonIndex(StructuralRegion::CDR2) <<
                ", CDR3: " << cdr_analyser.ShannonIndex(StructuralRegion::CDR3));
        INFO("Simpson index. CDR1: " << cdr_analyser.SimpsonIndex(StructuralRegion::CDR1) <<
             ", CDR2: " << cdr_analyser.SimpsonIndex(StructuralRegion::CDR2) <<
             ", CDR3: " << cdr_analyser.SimpsonIndex(StructuralRegion::CDR3));
        INFO("Clonal Shannon index: " << cdr_analyser.ClonalShannonIndex());
        INFO("Clonal Simpson index: " << cdr_analyser.ClonalSimpsonIndex());
        INFO("Diversity Analyzer ends");
    }

    void CDRLabelerLaunch::CheckInputParams() {
        const auto columns = ReportColumns::ColumnSet<DivanReportEvalContext>::ChooseColumns(
                config_.output_params.feature_report_params.preset, config_.output_params.feature_report_params.columns
        );
    }
}