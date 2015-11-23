#include "left_j_tail_aligner.hpp"

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>

using namespace std;
using namespace seqan;

IgGeneAlignmentPtr LeftJTailAligner::ComputeAlignment(IgGeneAlignmentPositions alignment_positions) {
    Align<Dna5String> align;
    resize(rows(align), 2);
    if(alignment_positions.alignment.query_pos.first == 1)
        return IgGeneAlignmentPtr(new IgGeneAlignment(alignment_positions, align));
    size_t tail_length = alignment_positions.alignment.subject_pos.first - 1;
    cout << "Tail length: " << tail_length << endl;
    auto read_segment = prefix(suffix(alignment_positions.read->seq, alignment_positions.alignment.query_pos.first -
            tail_length - left_shift_ - 1), tail_length + left_shift_);
    assignSource(row(align, 0), read_segment);
    assignSource(row(align, 1), prefix(alignment_positions.ig_gene->seq(), tail_length));
    globalAlignment(align, Score<int, Simple>(0, -1, -1));
    return IgGeneAlignmentPtr(new IgGeneAlignment(alignment_positions, align));
}