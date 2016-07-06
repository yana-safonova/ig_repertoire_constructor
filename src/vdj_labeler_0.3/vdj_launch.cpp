#include "vdj_launch.hpp"

#include <read_archive.hpp>
#include <germline_db_generator.hpp>
#include <vj_parallel_processor.hpp>
#include <block_alignment/block_alignment_converter.hpp>
#include <vdj_alignments/vdj_hits.hpp>
#include <vdj_alignments/vdj_hits_storage.hpp>
#include <alignment_utils/alignment_positions.hpp>
#include <vdj_alignments/aligners/simple_d_aligner.hpp>

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
    INFO("Generation of DB for join segments...");
    INFO("Alignment against VJ germline segments");
    vj_finder::VJParallelProcessor processor(read_archive, config_.vj_finder_config.algorithm_params,
                                             v_db, j_db,
                                             config_.run_params.threads_count);
    vj_finder::VJAlignmentInfo alignment_info = processor.Process();
    INFO(alignment_info.NumVJHits() << " reads were aligned; " << alignment_info.NumFilteredReads() <<
        " reads were filtered out");

    // INFO(alignment_info.AlignmentRecords()[6].Read());
    // INFO(alignment_info.AlignmentRecords()[6].VHits()[0].ImmuneGene());
    // INFO(alignment_info.AlignmentRecords()[0].VHits()[0].Read());
    // INFO(alignment_info.AlignmentRecords()[0].JHits().size());
    // auto vdj_hits = VDJHits(alignment_info.AlignmentRecords()[6]);
    // INFO((*(vdj_hits.VHits().cbegin()))->Alignment());

    auto vdj_storage = VDJHitsStorage(alignment_info);
    // INFO(*(vdj_storage[0]->Read()));

    alignment_utils::AlignmentPositions alignment_positions(std::make_pair<size_t, size_t>(100, read_archive[0].length() - 1),
                                                            std::make_pair<size_t, size_t>(0, 10));

    alignment_utils::ImmuneGeneAlignmentPositions immune_alignment_positions(alignment_positions,
                                                                             d_db[0],
                                                                             read_archive[0]);

    INFO(SimpleDAligner().ComputeAlignment(immune_alignment_positions)->Alignment());
}

} // End namespace vdj_labeler