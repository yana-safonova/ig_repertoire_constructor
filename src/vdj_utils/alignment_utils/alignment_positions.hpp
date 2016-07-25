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

    AlignmentPositions(const AlignmentPositions&) = default;
    AlignmentPositions(AlignmentPositions&&) = default;
    AlignmentPositions& operator=(const AlignmentPositions&) = default;
    AlignmentPositions& operator=(AlignmentPositions&&) = default;

    size_t QueryAlignmentLength() const { return query_pos.second - query_pos.first + 1; }

    size_t SubjectAlignmentLength() const { return subject_pos.second - subject_pos.first + 1; }
};

std::ostream &operator<<(std::ostream &out, const AlignmentPositions &obj);

// ------------------------------------------------------------
// Immune Gene Alignment Position
// ------------------------------------------------------------
class ImmuneGeneAlignmentPositions {
private:
    AlignmentPositions alignment_positions_;
    const germline_utils::ImmuneGene* immune_gene_ptr_;
    const core::Read* read_ptr_;

public:
    ImmuneGeneAlignmentPositions(AlignmentPositions alignment_positions,
                                 const germline_utils::ImmuneGene* immune_gene_ptr = nullptr,
                                 const core::Read* read_ptr = nullptr) :
            alignment_positions_(std::move(alignment_positions)),
            immune_gene_ptr_(immune_gene_ptr),
            read_ptr_(read_ptr)
    { }

    ImmuneGeneAlignmentPositions(AlignmentPositions alignment_positions,
                                 const germline_utils::ImmuneGene &immune_gene,
                                 const core::Read &read) :
        ImmuneGeneAlignmentPositions(std::move(alignment_positions), &immune_gene, &read)
    { }

    ImmuneGeneAlignmentPositions(const ImmuneGeneAlignmentPositions&) = default;
    ImmuneGeneAlignmentPositions(ImmuneGeneAlignmentPositions&&) = default;
    ImmuneGeneAlignmentPositions& operator=(const ImmuneGeneAlignmentPositions&) = default;
    ImmuneGeneAlignmentPositions& operator=(ImmuneGeneAlignmentPositions&&) = default;

    size_t GeneStartPos() const { return alignment_positions_.subject_pos.first; }
    size_t GeneEndPos  () const { return alignment_positions_.subject_pos.second; }
    size_t ReadStartPos() const { return alignment_positions_.query_pos.first; }
    size_t ReadEndPos  () const { return alignment_positions_.query_pos.second; }

    bool IsEmpty() const { return ReadStartPos() > ReadEndPos(); }

    size_t ReadAlignmentLength() const { return alignment_positions_.QueryAlignmentLength(); }

    size_t GeneAlignmentLength() const { return alignment_positions_.SubjectAlignmentLength(); }

    const germline_utils::ImmuneGene* GenePtr() const { return immune_gene_ptr_; }
    const germline_utils::ImmuneGene& Gene() const {
        assert(immune_gene_ptr_ != nullptr);
        return *immune_gene_ptr_;
    }

    const core::Read* ReadPtr() const { return read_ptr_; }
    const core::Read& Read() const {
        assert(read_ptr_ != nullptr);
        return *read_ptr_;
    }

    const AlignmentPositions &AlignmentPositions() const { return alignment_positions_; }
};

std::ostream& operator<<(std::ostream& out, const ImmuneGeneAlignmentPositions& obj);
} // End namespace alignment_utils