#include <logger/logger.hpp>
#include <block_alignment/pairwise_block_alignment.hpp>
#include <cdr3_hamming_graph_connected_components_processors/base_cdr3_hg_cc_processor.hpp>
#include "parent_read_reconstructor.hpp"

namespace antevolo {

    std::tuple<core::Read,
            seqan::Align<seqan::Dna5String>,
            seqan::Align<seqan::Dna5String>> ParentReadReconstructor::ReconstructParentRead(
           const annotation_utils::AnnotatedClone& clone1,
           const annotation_utils::AnnotatedClone& clone2,
           size_t id,
           size_t gene_cdr3_start_pos,
           size_t gene_cdr3_end_pos) {


        VERIFY_MSG(Base_CDR3_HG_CC_Processor::CheckClonesConsistencyForReconstruction(clone1, clone2),
                   "this method can be used only if " <<
                   "Base_CDR3_HG_CC_Processor::CheckClonesConsistencyForReconstruction returns true");
        // !!
        // can be used only if
        // Base_CDR3_HG_CC_Processor::CheckClonesConsistencyForReconstruction == true

        // read
        std::stringstream s;
        s << std::string("fake_read_") << id << ("___size___1");
        std::string read_name = s.str();
        std::string res_string;
        seqan::Dna5String res_dna5string;

        auto read1 = clone1.Read().seq;
        auto read2 = clone2.Read().seq;

        const auto& v_gene = clone1.VGene();
        const auto& j_gene = clone1.JGene();

        const TRow& gene_v_alignment1 = seqan::row(clone1.VAlignment().Alignment(), 0);
        const TRow& gene_v_alignment2 = seqan::row(clone2.VAlignment().Alignment(), 0);
        const TRow& v_alignment1 = seqan::row(clone1.VAlignment().Alignment(), 1);
        const TRow& v_alignment2 = seqan::row(clone2.VAlignment().Alignment(), 1);
        size_t v_gene_cdr3_start_view_pos1 = seqan::toViewPosition(
                gene_v_alignment1,
                gene_cdr3_start_pos);
        size_t v_gene_cdr3_start_view_pos2 = seqan::toViewPosition(
                gene_v_alignment2,
                gene_cdr3_start_pos);

        // since insertion blocks are equal,
        // gene alignments must be identical before CDR3 start pos,
        // check for equality
        {
            bool eq = true;
            VERIFY(v_gene_cdr3_start_view_pos1 == v_gene_cdr3_start_view_pos2);
            for (size_t i = 0; i < v_gene_cdr3_start_view_pos1; ++i) {
                if (gene_v_alignment1[i] !=
                    gene_v_alignment2[i]) {
                    eq = false;
                    break;
                }
            }
            VERIFY(eq);
        }
        size_t alignment_length = length(v_alignment1);
        // [0, cdr3_start_pos)
        auto read_v_gaps = TraverseAlignments(v_alignment1,
                                              v_alignment2,
                                              gene_v_alignment1,
                                              0,
                                              0,
                                              v_gene_cdr3_start_view_pos1,
                                              v_gene_cdr3_start_view_pos2,
                                              res_string);
        // [cdr3_start_pos, v_gene_end_pos]
//        for (size_t i = toViewPosition(v_alignment1, clone1.CDR3Range().start_pos); i < alignment_length; ++i) {
//        for (size_t i = toViewPosition(v_alignment1, clone1.CDR3Range().start_pos - 1) + 1; i < alignment_length; ++i) {
        for (size_t i = res_string.length() + read_v_gaps.size(); i < alignment_length; ++i) {

            if (v_alignment1[i] == '-') {
                read_v_gaps.push_back(i);
            }
            else {
                res_string.push_back(v_alignment1[i]);
            }
        }
        // (v_gene_end_pos, j_gene_start_pos)
        // -- the CDR3 part outside V and J
        size_t v_clipping_end_position = res_string.length() + read_v_gaps.size();
        for (size_t i = clone1.VAlignment().EndQueryPosition() + 1;
             i < clone1.JAlignment().StartQueryPosition();
             ++i) {
            res_string.push_back(read1[i]);
        }
        size_t j_clipping_begin_position;
        if (clone1.VAlignment().EndQueryPosition() < clone1.JAlignment().StartQueryPosition()) {
            j_clipping_begin_position = res_string.length();
        }
        else {
            j_clipping_begin_position = res_string.length() - 1 -
                    clone1.VAlignment().EndQueryPosition() + clone1.JAlignment().StartQueryPosition();
        }
        
        const TRow& gene_j_alignment1 = seqan::row(clone1.JAlignment().Alignment(), 0);
        const TRow& gene_j_alignment2 = seqan::row(clone2.JAlignment().Alignment(), 0);
        const TRow& j_alignment1 = seqan::row(clone1.JAlignment().Alignment(), 1);
        const TRow& j_alignment2 = seqan::row(clone2.JAlignment().Alignment(), 1);
        size_t j_gene_cdr3_end_view_pos1 = seqan::toViewPosition(
                gene_j_alignment1,
                gene_cdr3_end_pos);
        size_t j_gene_cdr3_end_view_pos2 = seqan::toViewPosition(
                gene_j_alignment2,
                gene_cdr3_end_pos);
        // again check for gene alignments' equality
        {
            bool eq = true;
            size_t al_length = length(gene_j_alignment1);
            VERIFY(al_length - j_gene_cdr3_end_view_pos1 == length(gene_j_alignment2) - j_gene_cdr3_end_view_pos2);
            for (size_t i = 1; i < al_length - j_gene_cdr3_end_view_pos1; ++i) {
                if (gene_j_alignment1[i + j_gene_cdr3_end_view_pos1] !=
                    gene_j_alignment2[i + j_gene_cdr3_end_view_pos2]) {
                    eq = false;
                    break;
                }
            }
            VERIFY(eq);
        }
        std::vector<size_t> read_j_gaps;
        // [j_gene_start_pos, cdr3_end_pos]
        for (size_t i = 0;
             i <= j_gene_cdr3_end_view_pos1;
             ++i) {
            if (j_alignment1[i] == '-') {
                read_j_gaps.push_back(i);
            }
            else {
                if (toSourcePosition(j_alignment1, i) > clone1.VAlignment().EndQueryPosition()) {
                    // there can be nothing between V and J
                    res_string.push_back(j_alignment1[i]);
                }
            }
        }
        // (cdr3_end_pos, read_end_pos]
        auto pre_j_gaps = TraverseAlignments(j_alignment1,
                                             j_alignment2,
                                             gene_j_alignment1,
                                             j_gene_cdr3_end_view_pos1 + 1,
                                             j_gene_cdr3_end_view_pos2 + 1,
                                             length(j_alignment1),
                                             length(j_alignment2),
                                             res_string);
        for (size_t g : pre_j_gaps) {
            read_j_gaps.push_back(g);
        }

        for (auto c : res_string) {
            seqan::append(res_dna5string, c);
        }
        core::Read res_read(read_name, res_dna5string, id);

        using namespace seqan;
        Align<Dna5String, ArrayGaps> seqan_v_alignment;
        resize(rows(seqan_v_alignment), 2);
        assignSource(row(seqan_v_alignment, 0), v_gene.seq());
        assignSource(row(seqan_v_alignment, 1), res_read.seq);
        TRow& v_gene_row = row(seqan_v_alignment, 0);
        TRow& v_read_row = row(seqan_v_alignment, 1);
        v_gene_row = gene_v_alignment1;
        for (size_t i : read_v_gaps) {
            insertGap(v_read_row, i);
        }

        setClippedEndPosition(v_read_row, v_clipping_end_position);

        Align<Dna5String, ArrayGaps> seqan_j_alignment;
        resize(rows(seqan_j_alignment), 2);
        assignSource(row(seqan_j_alignment, 0), j_gene.seq());
        assignSource(row(seqan_j_alignment, 1), res_read.seq);
        TRow& j_gene_row = row(seqan_j_alignment, 0);
        TRow& j_read_row = row(seqan_j_alignment, 1);
        j_gene_row = gene_j_alignment1;

        setClippedBeginPosition(j_read_row, j_clipping_begin_position);
        for (size_t i : read_j_gaps) {
            insertGap(j_read_row, i);
        }


        if (length(v_gene_row) != length(v_read_row)) {
            INFO(" V: " << length(v_gene_row) << "\t" << length(v_read_row));
        }
        if (length(j_gene_row) != length(j_read_row)) {
            INFO(" J: " << length(j_gene_row) << "\t" << length(j_read_row));
        }
        VERIFY(length(v_gene_row) == length(v_read_row) && length(j_gene_row) == length(j_read_row));

        std::tuple<core::Read, Align<Dna5String, ArrayGaps>, Align<Dna5String, ArrayGaps>> res(
                res_read,
                seqan_v_alignment,
                seqan_j_alignment);
        return res;
    }




    /* !!
     *  this function can only be applied to reads
     *  with equal insertion blocks
    */
    std::vector<size_t> ParentReadReconstructor::TraverseAlignments(
            const TRow& alignment1,
            const TRow& alignment2,
            const TRow& gene_alignment,
            size_t view1_start_pos,
            size_t view2_start_pos,
            size_t view1_end_pos,
            size_t view2_end_pos,
            std::string &res_string) {

        using namespace seqan;
        std::vector<size_t> read_gaps;

        size_t align1_start_pos = view1_start_pos;
        size_t align2_start_pos = view2_start_pos;
        size_t align1_end_pos = view1_end_pos;
        size_t align2_end_pos = view2_end_pos;

        VERIFY_MSG(align1_end_pos  - align1_start_pos == align2_end_pos - align2_start_pos,
                   align1_start_pos << "\t" << align2_start_pos
                                    << "\t" << align1_end_pos << "\t" << align2_end_pos);

        size_t i = align1_start_pos;
        size_t j = align2_start_pos;
        while (i <  align1_end_pos &&
               j <  align2_end_pos) {
            if (alignment1[i] == alignment2[j]) {
                if (alignment1[i] == '-') {
                    read_gaps.push_back(i);
                }
                else {
                    res_string.push_back(alignment1[i]);
                }
            }
            else {
                VERIFY(gene_alignment[i] != '-')
                res_string.push_back(gene_alignment[i]);
            }
            ++i;
            ++j;

        }
        return read_gaps;
    }
}