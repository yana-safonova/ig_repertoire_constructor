#include "vdj_launch.hpp"

#include <read_archive.hpp>
#include <germline_db_generator.hpp>
#include <vj_parallel_processor.hpp>
#include <germline_utils/germline_databases/chain_database.hpp>
#include <vdj_alignments/vdj_hits/vdj_hits.hpp>
#include <vdj_alignments/vdj_hits/single_immune_gene_segment_hits.hpp>
#include <vdj_alignments/vdj_hits/d_gene_segment_hit.hpp>
#include <vdj_alignments/vdj_hits/vdj_hits_storage.hpp>
// #include <alignment_utils/alignment_positions.hpp>
#include <vdj_alignments/aligners/simple_d_aligner.hpp>
#include <vdj_alignments/hits_calculator/alignment_quality_checkers/match_threshold_alignment_quality_checker.hpp>
#include <vdj_alignments/hits_calculator/d_alignment_positions_checkers/info_based_d_alignment_position_checker.hpp>
#include <vdj_alignments/hits_calculator/d_alignment_positions_calculator/custom_d_alignment_positions_calculator.hpp>
#include <vdj_alignments/hits_calculator/d_hits_calculator/info_based_d_hits_calculator.hpp>
// #include <vdj_alignments/hits_calculator/alignment_quality_checkers/threshold_alignment_estimator.hpp>

namespace vdj_labeler {

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

    SimpleDAligner d_aligner;
    MatchThresholdAlignmentQualityChecker quality_checker(7);
    // ThresholdAlignmentQualityChecker quality_checker(1);
    InfoBasedDAlignmentPositionChecker position_checker(config_.d_align_quality_params);
    CustomDAlignmentPositionsCalculator positions_calculator;
    InfoBasedDHitsCalculator calculator(
        d_db.GetConstDbByGeneType(ImmuneGeneType(ChainType("IGH"), SegmentType::DiversitySegment)),
        d_aligner, quality_checker, position_checker, positions_calculator,
        config_.d_align_quality_params);

    VDJHitsStorage vdj_storage (alignment_info);
    VDJHitsStorage vdj_storage2(alignment_info, calculator);

    std::vector<size_t> n_d_hits;
    n_d_hits.reserve(vdj_storage2.size());
    for (const auto& vdj_hit : vdj_storage2) {
        if (vdj_hit.DHits().size() == 0) {
            n_d_hits.push_back(0);
        }
        for (const auto& d_hits : vdj_hit.DHits()) {
            if (d_hits.size() >= 3) {
                INFO(vdj_hit.Read().seq);
                for (const auto& d_hit : d_hits) {
                    INFO(d_hit.Subject().name());
                }
            }
            n_d_hits.push_back(d_hits.size());
        }
    }
    std::vector<size_t> stat_n_d_hits(10);
    for (auto n : n_d_hits) {
        stat_n_d_hits[n]++;
    }

    for (auto freq : stat_n_d_hits) {
        INFO(freq);
    }

    // //TestRecombinationCalculator(read_archive, vdj_storage2);


    // size_t max_cleavage = 20;
    // size_t max_palindrome = 7;
    // LeftEventSHMsCalculator left_shms_calculator;
    // RightEventSHMsCalculator right_shms_calculator;
    // VersatileGeneSHMsCalculator shms_calculator(left_shms_calculator, right_shms_calculator);
    // VRecombinationEventGenerator v_generator(shms_calculator, max_cleavage, max_palindrome);
    // DRecombinationEventGenerator d_generator(shms_calculator, max_cleavage, max_palindrome);
    // JRecombinationEventGenerator j_generator(shms_calculator, max_cleavage, max_palindrome);
    // VersatileInsertionGenerator insertion_generator;
    // CustomHeavyChainRecombinationGenerator recombination_generator(v_generator,
    //                                                                d_generator,
    //                                                                j_generator,
    //                                                                insertion_generator,
    //                                                                insertion_generator);
    // HcRecombinationEstimator recombination_estimator;
    // INFO("Generator of VDJ recombinations starts");
    // for(auto it = vdj_storage2.cbegin(); it != vdj_storage2.cend(); it++) {
    //     INFO(it->Read().name);
    //     INFO(it->VHits().size());
    //     INFO(it->DHits().size());
    //     INFO(it->JHits().size());
    //     INFO("");
    //     auto recombination_storage = recombination_generator.ComputeRecombinations(*it);
    //     // recombination_estimator.Update(recombination_storage);
    // }
    // INFO("Generator of VDJ recombinations ends");
    // recombination_estimator.OutputRecombinationNumber();
    // recombination_estimator.OutputSHMsDistribution();
    // recombination_estimator.OutputRecombinationEvents();
}

} // End namespace vdj_labeler