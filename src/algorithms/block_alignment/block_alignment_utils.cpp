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
        int finishing_gap = (gene_len - path.last().needle_pos) - (read_len - path.last().read_pos);
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
                const auto &gene_edge = infix(gene, path[i].needle_pos + path[i].length, path[i + 1].needle_pos);

                // INFO("EDGE: " << read_edge << " - " << gene_edge);

                if (length(read_edge) < length(gene_edge)) {
                    insertGaps(row_read,
                               path[i].read_pos + path[i].length + find_simple_gap(read_edge, gene_edge),
                               length(gene_edge) - length(read_edge));
                } else if (length(gene_edge) < length(read_edge)) {
                    insertGaps(row_gene,
                               path[i].needle_pos + path[i].length + find_simple_gap(gene_edge, read_edge),
                               length(read_edge) - length(gene_edge));
                } else {
                    // Don't add gaps
                }
            }
        }

        // Add starting gaps (reverse order since insert_gaps works with VIEW position!)
        int starting_gap = path.first().needle_pos - path.first().read_pos;
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

    // todo: add description
    std::vector<Match> combine_sequential_kmer_matches(std::vector<KmerMatch> &matches,
                                                       size_t K) {
        std::sort(matches.begin(), matches.end(), KmerMatch::less_shift);

        std::vector<Match> res;
        res.reserve(matches.size()); // TODO Is it really necessary?

        if (matches.size() == 0) {
            return res;
        }

        Match cur = { matches[0].needle_pos, matches[0].read_pos, K }; // start first match
        for (size_t i = 1; i < matches.size(); ++i) {
            if (matches[i].needle_pos == matches[i-1].needle_pos + 1 &&
                    matches[i].read_pos == matches[i-1].read_pos + 1) { // extend current match
                cur.length += 1;
            } else { // save match and start new one
                res.push_back(cur);
                cur = { matches[i].needle_pos, matches[i].read_pos, K };
            }
        }
        res.push_back(cur); // save last match
        return res;
    }
}