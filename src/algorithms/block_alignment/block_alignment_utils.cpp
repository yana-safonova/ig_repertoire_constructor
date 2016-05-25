#include <verify.hpp>
#include "block_alignment_utils.hpp"

namespace algorithms {
    // function converts alignment path into Alignment object
    // TODO: make it consistent with other code
    TAlign path2seqanAlignment(const AlignmentPath &path, const seqan::Dna5String &read, const seqan::Dna5String &gene) {
        VERIFY(path.check_overlaps());
        VERIFY(!path.empty());

        using namespace seqan;
        using TRow =  seqan::Row<TAlign>::Type; // gapped sequence type

        TAlign align;
        resize(rows(align), 2);
        assignSource(row(align, 0), gene);
        assignSource(row(align, 1), read);

        TRow & row_gene = row(align, 0);
        TRow & row_read = row(align, 1);

        // Clip head
        size_t read_len = length(read);
        size_t gene_len = length(gene);

        // Add finishing gaps (reverse order since insert_gaps works with VIEW position!)
        // TODO Use SeqAn clipping if possible
        int finishing_gap = (int(gene_len) - path.last().subject_pos) - (int(read_len) - path.last().read_pos);
        if (finishing_gap > 0) {
            insertGaps(row_read, read_len, finishing_gap);
        } else if (finishing_gap < 0) {
            // insertGaps(row_gene, gene_len, -finishing_gap);
            // Use clipping:
            setEndPosition(row_read, read_len + finishing_gap);
        } else {
            // Do nothing
        }

        // Add edge gaps if needed. Reverse order!
        if (path.size() > 1) {
            for (size_t i = path.size() - 2; i + 1 > 0; --i) {
                const auto &read_edge = infix(read, path[i].read_pos + path[i].length, path[i + 1].read_pos);
                const auto &gene_edge = infix(gene, path[i].subject_pos + path[i].length, path[i + 1].subject_pos);

                // INFO("EDGE: " << read_edge << " - " << gene_edge);

                if (length(read_edge) < length(gene_edge)) {
                    insertGaps(row_read,
                               path[i].read_pos + path[i].length + find_simple_gap(read_edge, gene_edge),
                               length(gene_edge) - length(read_edge));
                } else if (length(gene_edge) < length(read_edge)) {
                    insertGaps(row_gene,
                               path[i].subject_pos + path[i].length + find_simple_gap(gene_edge, read_edge),
                               length(read_edge) - length(gene_edge));
                } else {
                    // Don't add gaps
                }
            }
        }

        // Add starting gaps (reverse order since insert_gaps works with VIEW position!)
        int starting_gap = path.first().subject_pos - path.first().read_pos;
        if (starting_gap > 0) {
            insertGaps(row_read, 0, starting_gap);
        } else if (starting_gap < 0) {
            // insertGaps(row_gene, 0, -starting_gap);
            // Use clipping:
            setBeginPosition(row_read, -starting_gap);
        } else {
            // Do nothing
        }

        return align;
    }
}