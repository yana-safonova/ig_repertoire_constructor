#include "vdj_launch.hpp"

#include <read_archive.hpp>
#include <germline_db_generator.hpp>
#include <vj_parallel_processor.hpp>
#include <block_alignment/block_alignment_converter.hpp>
#include <alignment_utils/pairwise_alignment.hpp>
#include <block_alignment/block_alignment_converter.hpp>

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
    CustomGeneDatabase j_db(SegmentType::JoinSegment);

    v_db.AddDatabase(ImmuneGeneType(ChainType("IGH"), SegmentType::VariableSegment), v_germline_genes_fname);
    j_db.AddDatabase(ImmuneGeneType(ChainType("IGH"), SegmentType::JoinSegment), j_germline_genes_fname);
    INFO("Generation of DB for join segments...");
    INFO("Alignment against VJ germline segments");
    vj_finder::VJParallelProcessor processor(read_archive, config_.vj_finder_config.algorithm_params,
                                             v_db, j_db,
                                             config_.run_params.threads_count);
    vj_finder::VJAlignmentInfo alignment_info = processor.Process();
    INFO(alignment_info.NumVJHits() << " reads were aligned; " << alignment_info.NumFilteredReads() <<
        " reads were filtered out");

    // algorithms::BlockAlignmentConverter block_alignment_converter;
    // for (auto& alignment_record : alignment_info.AlignmentRecords()) {
    //     block_alignment_converter.ConvertToAlignment(alignment_record.)
    // }
}

}