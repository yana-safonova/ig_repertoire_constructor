#include "vdj_launch.hpp"

#include <read_archive.hpp>
#include <germline_db_generator.hpp>
#include <vj_parallel_processor.hpp>
#include <block_alignment/block_alignment_converter.hpp>
#include <vdj_alignments/vdj_hits.hpp>
#include <vdj_alignments/vdj_hits_storage.hpp>
#include <alignment_utils/alignment_positions.hpp>
#include <vdj_alignments/aligners/simple_d_aligner.hpp>
#include <model/recombination_model.hpp>
#include <germline_utils/chain_type.hpp>
#include <recombination_calculator/hc_model_based_recombination_calculator.hpp>
#include <recombination_utils/recombination_storage.hpp>
#include <vdj_alignments/hits_calculator/d_alignment_positions_checkers/info_based_d_alignment_position_checker.hpp>
#include <vdj_alignments/hits_calculator/d_alignment_positions_calculator/custom_d_alignment_positions_calculator.hpp>
#include <vdj_alignments/hits_calculator/d_hits_calculator/info_based_d_hits_calculator.hpp>
#include <vdj_alignments/hits_calculator/alignment_quality_checkers/match_threshold_alignment_quality_checker.hpp>
#include <vdj_alignments/hits_calculator/alignment_quality_checkers/threshold_alignment_estimator.hpp>
#include <recombination_generation/gene_events_generators/shms_calculators/left_event_shms_calculator.hpp>
#include <recombination_generation/gene_events_generators/shms_calculators/right_event_shms_calculator.hpp>
#include <recombination_generation/gene_events_generators/shms_calculators/versatile_shms_calculator.hpp>
#include <recombination_generation/gene_events_generators/v_recombination_event_generator.hpp>
#include <recombination_generation/gene_events_generators/d_recombination_event_generator.hpp>
#include <recombination_generation/gene_events_generators/j_recombination_event_generator.hpp>
#include <recombination_generation/insertion_events_generators/versatile_insertion_event_generator.hpp>
#include <recombination_generation/custom_hc_recombination_generator.hpp>
#include <recombination_estimators/hc_recombination_estimator.hpp>


namespace vdj_labeler {

// void VDJLabelerLaunch::TestRecombinationCalculator(const core::ReadArchive& reads_archive,
//                                                    const VDJHitsStorage &hits_storage)
// {
//     size_t read_index = 3;
//     core::ReadPtr read_3 = std::make_shared<core::Read>(reads_archive[read_index]);
//     VDJHitsPtr hits_3 = hits_storage[read_index];
//     INFO("Read 3. #V: " << hits_3->VHitsNumber() <<
//         ", #D: " << hits_3->DHitsNumber() <<
//         ", #J: " << hits_3->JHitsNumber());
//
//     auto v_alignment = hits_3->GetAlignmentByIndex(germline_utils::SegmentType::VariableSegment, 0);
//     recombination_utils::CleavedIgGeneAlignment v_event_0(v_alignment, 0, 0, 0, 0);
//     recombination_utils::CleavedIgGeneAlignment v_event_1(v_alignment, 0, -1, 0, 0);
//     recombination_utils::CleavedIgGeneAlignment v_event_2(v_alignment, 0, -2, 0, 0);
//     recombination_utils::CleavedIgGeneAlignment v_event_3(v_alignment, 0, -3, 0, 1);
//
//     auto d_alignment = hits_3->GetAlignmentByIndex(germline_utils::SegmentType::DiversitySegment, 0);
//     recombination_utils::CleavedIgGeneAlignment d_event_0(d_alignment, 1, 8, 0, 0);
//
//     auto j_alignment = hits_3->GetAlignmentByIndex(germline_utils::SegmentType::JoinSegment, 0);
//     recombination_utils::CleavedIgGeneAlignment j_event_0(j_alignment, 0, 0, 1, 0);
//     recombination_utils::CleavedIgGeneAlignment j_event_1(j_alignment, 1, 0, 0, 0);
//
//     recombination_utils::NongenomicInsertion vd_insertion_0(425, 441);
//     recombination_utils::NongenomicInsertion vd_insertion_1(426, 441);
//     recombination_utils::NongenomicInsertion vd_insertion_2(427, 441);
//     recombination_utils::NongenomicInsertion vd_insertion_3(428, 441);
//
//     recombination_utils::NongenomicInsertion dj_insertion_0(453, 452);
//     recombination_utils::NongenomicInsertion dj_insertion_1(453, 453);
//
//     recombination_utils::RecombinationStorage<recombination_utils::HCRecombination> recombination_storage(read_3);
//     recombination_storage.AddRecombination(recombination_utils::HCRecombination(read_3, v_event_0, d_event_0, j_event_0,
//                                                            vd_insertion_0, dj_insertion_0));
//     recombination_storage.AddRecombination(recombination_utils::HCRecombination(read_3, v_event_1, d_event_0, j_event_0,
//                                                            vd_insertion_1, dj_insertion_0));
//     recombination_storage.AddRecombination(recombination_utils::HCRecombination(read_3, v_event_2, d_event_0, j_event_0,
//                                                            vd_insertion_2, dj_insertion_0));
//     recombination_storage.AddRecombination(recombination_utils::HCRecombination(read_3, v_event_3, d_event_0, j_event_0,
//                                                            vd_insertion_3, dj_insertion_0));
//     recombination_storage.AddRecombination(recombination_utils::HCRecombination(read_3, v_event_0, d_event_0, j_event_1,
//                                                            vd_insertion_0, dj_insertion_1));
//     recombination_storage.AddRecombination(recombination_utils::HCRecombination(read_3, v_event_1, d_event_0, j_event_1,
//                                                            vd_insertion_1, dj_insertion_1));
//     recombination_storage.AddRecombination(recombination_utils::HCRecombination(read_3, v_event_2, d_event_0, j_event_1,
//                                                            vd_insertion_2, dj_insertion_1));
//     recombination_storage.AddRecombination(recombination_utils::HCRecombination(read_3, v_event_3, d_event_0, j_event_1,
//                                                            vd_insertion_3, dj_insertion_1));
//     INFO(recombination_storage.size() << " recombinaions were generated");
// }

void VDJLabelerLaunch::Launch() {
    INFO("VDJ labeler starts");
    std::string input_filename = config_.io_params.input_params.input_sequences;
    std::string v_germline_genes_fname = config_.io_params.input_params.germlines.igh_genes.variable_genes;
    std::string d_germline_genes_fname = config_.io_params.input_params.germlines.igh_genes.diversity_genes;
    std::string j_germline_genes_fname = config_.io_params.input_params.germlines.igh_genes.join_genes;

    core::ReadArchive read_archive(input_filename);
    read_archive.FixSpacesInHeaders();

    using namespace germline_utils;
    CustomGeneDatabase v_db(SegmentType::VariableSegment);
    CustomGeneDatabase d_db(SegmentType::DiversitySegment);
    CustomGeneDatabase j_db(SegmentType::JoinSegment);

    v_db.AddDatabase(ImmuneGeneType(ChainType("IGH"), SegmentType::VariableSegment), v_germline_genes_fname);
    d_db.AddDatabase(ImmuneGeneType(ChainType("IGH"), SegmentType::DiversitySegment), d_germline_genes_fname);
    j_db.AddDatabase(ImmuneGeneType(ChainType("IGH"), SegmentType::JoinSegment), j_germline_genes_fname);

    germline_utils::ChainDatabase hc_db(germline_utils::ImmuneChainType::HeavyIgChain);
    hc_db.AddGenesFromFile(SegmentType::VariableSegment, v_germline_genes_fname);
    hc_db.AddGenesFromFile(SegmentType::DiversitySegment, d_germline_genes_fname);
    hc_db.AddGenesFromFile(SegmentType::JoinSegment, j_germline_genes_fname);

    INFO("Generation of DB for join segments...");
    INFO("Alignment against VJ germline segments");
    vj_finder::VJParallelProcessor processor(read_archive, config_.vj_finder_config.algorithm_params,
                                             v_db, j_db,
                                             config_.run_params.threads_count);
    vj_finder::VJAlignmentInfo alignment_info = processor.Process();
    INFO(alignment_info.NumVJHits() << " reads were aligned; " << alignment_info.NumFilteredReads() <<
        " reads were filtered out");

    // Andy: Blank model "tested" here
    // {
    //     std::ifstream in("src/vdj_labeler_0.3/test/blank_model.csv");
    //     // HCProbabilityRecombinationModel model(in, hc_db);
    //     // std::cout << model;
    //     IgGeneProbabilityModel model_V(in, hc_db.GetDb(germline_utils::SegmentType::VariableSegment));
    //     // std::cout << model_V;
    //     IgGeneProbabilityModel model_D(in, hc_db.GetDb(germline_utils::SegmentType::DiversitySegment));
    //     // std::cout << model_D;
    //     // IgGeneProbabilityModel model_J(in, hc_db.JoinGenes());
    //     // // cout << model_J;
    //     // NongenomicInsertionModel modelVD(in);
    //     // NongenomicInsertionModel modelDJ(in);
    //     // // cout << modelVD;
    //     // // cout << modelDJ;
    //     // PalindromeDeletionModel modelDelV(in, hc_db.VariableGenes());
    //     // // cout << modelDelV;
    //     // // cout << modelDelV.GetIgGeneDatabase() -> GetByIndex(0) -> name() << " ";
    //     // // cout << modelDelV.GetDeletionProbability(0, 0) << endl;
    //     // PalindromeDeletionModel modelDelJ(in, hc_db.JoinGenes());
    //     // // cout << modelDelJ.GetIgGeneDatabase() -> GetByIndex(1) -> name() << " ";
    //     // // cout << modelDelJ.GetDeletionProbability(1, -2) << endl;
    //     // PalindromeDeletionModel modelDelDL(in, hc_db.DiversityGenes());
    //     // PalindromeDeletionModel modelDelDR(in, hc_db.DiversityGenes());

    //     // HCModelBasedRecombinationCalculator recombination_calculator(model);
    // }

    SimpleDAligner d_aligner;
    MatchThresholdAlignmentQualityChecker quality_checker(5);
    // ThresholdAlignmentQualityChecker quality_checker(1);
    InfoBasedDAlignmentPositionChecker position_checker(config_.d_align_quality_params);
    CustomDAlignmentPositionsCalculator positions_calculator;
    InfoBasedDHitsCalculator calculator(
        d_db.GetConstDbByGeneType(ImmuneGeneType(ChainType("IGH"), SegmentType::DiversitySegment)),
        d_aligner, quality_checker, position_checker, positions_calculator);

    VDJHitsStorage vdj_storage (alignment_info);
    VDJHitsStorage vdj_storage2(alignment_info, calculator);

    //TestRecombinationCalculator(read_archive, vdj_storage2);


    size_t max_cleavage = 20;
    size_t max_palindrome = 7;
    LeftEventSHMsCalculator left_shms_calculator;
    RightEventSHMsCalculator right_shms_calculator;
    VersatileGeneSHMsCalculator shms_calculator(left_shms_calculator, right_shms_calculator);
    VRecombinationEventGenerator v_generator(shms_calculator, max_cleavage, max_palindrome);
    DRecombinationEventGenerator d_generator(shms_calculator, max_cleavage, max_palindrome);
    JRecombinationEventGenerator j_generator(shms_calculator, max_cleavage, max_palindrome);
    VersatileInsertionGenerator insertion_generator;
    CustomHeavyChainRecombinationGenerator recombination_generator(v_generator,
                                                                   d_generator,
                                                                   j_generator,
                                                                   insertion_generator,
                                                                   insertion_generator);
    HcRecombinationEstimator recombination_estimator;
    INFO("Generator of VDJ recombinations starts");
    for(auto it = vdj_storage2.cbegin(); it != vdj_storage2.cend(); it++) {
        INFO(it->Read().name);
        INFO(it->VHits().size());
        INFO(it->DHits().size());
        INFO(it->JHits().size());
        INFO("");
        auto recombination_storage = recombination_generator.ComputeRecombinations(*it);
        // recombination_estimator.Update(recombination_storage);
    }
    INFO("Generator of VDJ recombinations ends");
    recombination_estimator.OutputRecombinationNumber();
    recombination_estimator.OutputSHMsDistribution();
    recombination_estimator.OutputRecombinationEvents();
}

} // End namespace vdj_labeler