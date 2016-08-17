#pragma once

#include "pairwise_block_alignment.hpp"
#include <alignment_utils/pairwise_alignment.hpp>

namespace algorithms {
    template<typename SubjectTypename, typename QueryTypename>
    class BlockAlignmentConverter {

    protected:
        virtual seqan::Dna5String GetSubjectString(const SubjectTypename& subject) = 0;

        virtual seqan::Dna5String GetQueryString(const QueryTypename& query) = 0;

    public:
        alignment_utils::PairwiseAlignment<SubjectTypename, QueryTypename> ConvertToAlignment(
                const SubjectTypename &subject,
                const QueryTypename &query,
                PairwiseBlockAlignment block_alignment) {
            VERIFY_MSG(block_alignment.path.check_overlaps(), "Overlaps are not correct");
            VERIFY_MSG(!block_alignment.path.empty(), "Alignment path is empty");

            auto subject_string = GetSubjectString(subject); // in case of gene-read alignment, subject is gene
            auto query_string = GetQueryString(query); // in case of gene-read alignment, query is read

            using namespace seqan;
            using TRow =  seqan::Row<seqan::Align<Dna5String, seqan::ArrayGaps>>::Type; // gapped sequence type

            seqan::Align<Dna5String, seqan::ArrayGaps> align;
            resize(rows(align), 2);
            assignSource(row(align, 0), subject_string);
            assignSource(row(align, 1), query_string);

            TRow &row_gene = row(align, 0);
            TRow &row_read = row(align, 1);

            // Clip head
            size_t read_len = length(query_string);
            size_t gene_len = length(subject_string);

            // Add finishing gaps (reverse order since insert_gaps works with VIEW position!)
            // TODO Use SeqAn clipping if possible
            int finishing_gap = (int(gene_len) - block_alignment.path.last().subject_pos) -
                    (int(read_len) - block_alignment.path.last().read_pos);
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
            if (block_alignment.path.size() > 1) {
                for (size_t i = block_alignment.path.size() - 2; i + 1 > 0; --i) {
                    const auto &read_edge = infix(query_string,
                                                  block_alignment.path[i].read_pos + block_alignment.path[i].length,
                                                  block_alignment.path[i + 1].read_pos);
                    const auto &gene_edge = infix(subject_string,
                                                  block_alignment.path[i].subject_pos + block_alignment.path[i].length,
                                                  block_alignment.path[i + 1].subject_pos);

                    // INFO("EDGE: " << read_edge << " - " << gene_edge);

                    if (length(read_edge) < length(gene_edge)) {
                        insertGaps(row_read,
                                   block_alignment.path[i].read_pos + block_alignment.path[i].length +
                                           find_simple_gap(read_edge, gene_edge),
                                   length(gene_edge) - length(read_edge));
                    } else if (length(gene_edge) < length(read_edge)) {
                        insertGaps(row_gene,
                                   block_alignment.path[i].subject_pos + block_alignment.path[i].length +
                                           find_simple_gap(gene_edge, read_edge),
                                   length(read_edge) - length(gene_edge));
                    } else {
                        // Don't add gaps
                    }
                }
            }

            // Add starting gaps (reverse order since insert_gaps works with VIEW position!)
            int starting_gap = block_alignment.path.first().subject_pos - block_alignment.path.first().read_pos;
            if (starting_gap > 0) {
                insertGaps(row_read, 0, starting_gap);
            } else if (starting_gap < 0) {
                // insertGaps(row_gene, 0, -starting_gap);
                // Use clipping:
                setBeginPosition(row_read, -starting_gap);
            } else {
                // Do nothing
            }

            return alignment_utils::PairwiseAlignment<SubjectTypename, QueryTypename>(subject,
                                                                                      query,
                                                                                      align,
                                                                                      block_alignment.score);
        }
    };
}