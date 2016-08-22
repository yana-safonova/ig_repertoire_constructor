#pragma once

#include <read_archive.hpp>
#include "../germline_utils/germline_databases/immune_gene_database.hpp"

namespace alignment_utils {
    // ------------------------------------------------------------
    // Alignment Positions class
    // ------------------------------------------------------------
    // in our context, query is a reads, subject is a gene segment
    // both positions are inclusive
    struct AlignmentPositions {
        std::pair <size_t, size_t> query_pos;
        std::pair <size_t, size_t> subject_pos;

        AlignmentPositions() :
                query_pos(std::make_pair(size_t(-1), size_t(-1))),
                subject_pos(std::make_pair(size_t(-1), size_t(-1))) { }

        AlignmentPositions(std::pair <size_t, size_t> new_query_pos,
                           std::pair <size_t, size_t> new_subject_pos) :
                query_pos(new_query_pos),
                subject_pos(new_subject_pos) { }

        size_t QueryAlignmentLength() const { return query_pos.second - query_pos.first + 1; }

        size_t SubjectAlignmentLength() const { return subject_pos.second - subject_pos.first + 1; }
    };

    std::ostream &operator<<(std::ostream &out, const AlignmentPositions &obj);

    // ------------------------------------------------------------
    // Immune Gene Alignment Position
    // ------------------------------------------------------------
    class ImmuneGeneAlignmentPositions {
        AlignmentPositions alignment_positions_;
        const germline_utils::ImmuneGene& immune_gene_;
        const core::Read &read_;

//        ImmuneGeneAlignmentPositions() :
//                alignment_positions_(),
//                immune_gene_(),
//                read_(NULL) { }

        ImmuneGeneAlignmentPositions(AlignmentPositions alignment_positions,
                                     const germline_utils::ImmuneGene& immune_gene,
                                     const core::Read &read) :
                alignment_positions_(alignment_positions),
                immune_gene_(immune_gene),
                read_(read) { }

        size_t GeneStartPos() const { return alignment_positions_.subject_pos.first; }

        size_t GeneEndPos() const { return alignment_positions_.subject_pos.second; }

        size_t ReadStartPos() const { return alignment_positions_.query_pos.first; }

        size_t ReadEndPos() const { return alignment_positions_.query_pos.second; }

        bool IsEmpty() const { return ReadStartPos() > ReadEndPos(); }

        size_t ReadAlignmentLength() const { return alignment_positions_.QueryAlignmentLength(); }

        size_t GeneAlignmentLength() const { return alignment_positions_.SubjectAlignmentLength(); }

        const germline_utils::ImmuneGene& Gene() const { return immune_gene_; }

        const core::Read& Read() const { return read_; }
    };

    std::ostream& operator<<(std::ostream& out, const ImmuneGeneAlignmentPositions& obj);
}