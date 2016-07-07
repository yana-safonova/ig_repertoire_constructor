#include "logger/logger.hpp"
#include "simple_d_aligner.hpp"

using namespace seqan;

namespace vdj_labeler {

alignment_utils::ImmuneGeneReadAlignmentPtr SimpleDAligner::ComputeAlignment(
        alignment_utils::ImmuneGeneAlignmentPositions alignment_positions)
{
    TRACE("Computation of D alignment for positions: " << alignment_positions);
    Align<Dna5String> align;
    resize(rows(align), 2);

    if (alignment_positions.IsEmpty()) {
        TRACE("D alignment is empty");
        return std::make_shared<alignment_utils::ImmuneGeneReadAlignment>(alignment_positions.Gene(),
                                                                          alignment_positions.Read(),
                                                                          align);
    }

    auto read_segment = infixWithLength(alignment_positions.Read().seq,
                                        alignment_positions.ReadStartPos(),
                                        alignment_positions.ReadAlignmentLength());
    assignSource(row(align, 0), read_segment);
    assignSource(row(align, 1), alignment_positions.Gene().seq());
    auto& align_read = row(align, 0);
    TRACE("Read segment (" << length(read_segment) << "): " << read_segment);

    // Link: https://nar.oxfordjournals.org/content/early/2013/05/11/nar.gkt382.full.pdf
    // Citation:
    // SEARCH STRATEGY AND IMPLEMENTATION
    // ------ Default BLAST search parameters are used unless indicated otherwise;
    // ------ https://www.arabidopsis.org/Blast/BLASToptions.jsp
    // Identifying the V, D and J gene hits/:
    // ------ The default mismatch penalty is conservatively
    // ------ set to a relatively high value (-4) to minimize the
    // ------ chance of spurious matches.
    // Score::Score(match, mismatch, gap[, gapOpen]);
    double score = localAlignment(align, Score<int, Simple>(1, -4, -5, -2));
    // localAlignment(align, Score<int, Simple>(1, -4, 0, 0)); // For testing (allows to see many gaps)

    Align<Dna5String> align_orig(align);
    assignSource(row(align_orig, 0), alignment_positions.Read().seq);
    auto& align_origin_read = row(align_orig, 0);

    setClippedBeginPosition(align_origin_read, alignment_positions.ReadStartPos() + clippedBeginPosition(align_read));

    for(size_t read_pos = beginPosition(align_read);
        read_pos < length(align_read);
        ++read_pos)
    {
        if (isGap(align_read, read_pos))
            insertGap(align_origin_read, read_pos);
    }
    setEndPosition(align_origin_read, beginPosition(align_read) + length(align_read));

    return std::make_shared<alignment_utils::ImmuneGeneReadAlignment>(alignment_positions.Gene(),
                                                                      alignment_positions.Read(),
                                                                      align_orig,
                                                                      score);
}

} // End namespace vdj_labeler
