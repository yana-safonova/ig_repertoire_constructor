#include "logger/logger.hpp"
#include "info_based_d_hits_calculator.hpp"

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>

using namespace std;

namespace vdj_labeler {

alignment_utils::AlignmentPositions InfoBasedDHitsCalculator::ComputeDPositions(
        alignment_utils::ImmuneGeneAlignmentPositions v_positions,
        alignment_utils::ImmuneGeneAlignmentPositions j_positions)
{
    return alignment_utils::AlignmentPositions(make_pair(v_positions.AlignmentPositions().query_pos.second + 1,
                                                         j_positions.AlignmentPositions().query_pos.first - 1),
                                               make_pair(0, size_t(-1)));
}

alignment_utils::ImmuneGeneAlignmentPositions InfoBasedDHitsCalculator::CreateDAlignmentPositions(
        alignment_utils::AlignmentPositions d_alignment_positions,
        germline_utils::ImmuneGenePtr gene_ptr,
        core::ReadPtr read_ptr)
{
    d_alignment_positions.subject_pos.second = seqan::length(gene_ptr->seq()) - 1;
    return alignment_utils::ImmuneGeneAlignmentPositions(d_alignment_positions, *gene_ptr, *read_ptr);
}

// todo: refactor magic number!
bool InfoBasedDHitsCalculator::DAlignmentPositionsAreGood(alignment_utils::AlignmentPositions d_alignment_positions) {
    return d_alignment_positions.QueryAlignmentLength() >= 5;
}

// vdj_labeler::ImmuneGeneSegmentHitsPtr InfoBasedDHitsCalculator::ComputeHits(core::ReadPtr read_ptr) {
//     size_t read_index = read_archive_.GetIndexByReadName(read_ptr->name);
//     assert(read_index == read_ptr->id);
//
//     alignment_utils::ImmuneGeneAlignmentPositions v_alignment_positions =
//         vj_alignment_info_.GetVAlignmentByReadIndex(read_index);
//     alignment_utils::ImmuneGeneAlignmentPositions j_alignment_positions =
//         vj_alignment_info_.GetJAlignmentByReadIndex(read_index);
//     alignment_utils::AlignmentPositions d_positions = ComputeDPositions(v_alignment_positions, j_alignment_positions);
//
//     ImmuneGeneSegmentHitsPtr d_hits_ptr(new vdj_labeler::ImmuneGeneSegmentHits(
//                                                                         germline_utils::SegmentType::DiversitySegment,
//                                                                         read_ptr));
//     // if (!DAlignmentPositionsAreGood(d_positions)) {
//     //     TRACE("D positions are too short to generate alignment");
//     //     TRACE(d_positions);
//     //     // add single empty alignment and return hits storage
//     //     seqan::Align<seqan::Dna5String> align;
//     //     seqan::resize(seqan::rows(align), 2);
//     //     d_hits_ptr->AddHit(alignment_utils::ImmuneGeneReadAlignmentPtr(new alignment_utils::ImmuneGeneReadAlignment(
//     //         CreateDAlignmentPositions(d_positions,
//     //                                   make_shared<germline_utils::ImmuneGene>(*d_gene_database_.cbegin()),
//     //                                   read_ptr),
//     //         align, -1)));
//     //     return d_hits_ptr;
//     // }
//     // for (auto d_gene = d_gene_database_.cbegin(); d_gene != d_gene_database_.cend(); d_gene++) {
//     //     alignment_utils::ImmuneGeneAlignmentPositions d_alignment_pos = CreateDAlignmentPositions(d_positions,
//     //                                                                  make_shared<germline_utils::ImmuneGene>(*d_gene),
//     //                                                                  read_ptr);
//     //     auto d_alignment = d_gene_aligner_.ComputeAlignment(d_alignment_pos);
//     //     if (quality_checker_.AlignmentIsGood(d_alignment))
//     //         d_hits_ptr->AddHit(d_alignment);
//     // }
    // return d_hits_ptr;
// }

} // End namespace vdj_labeler