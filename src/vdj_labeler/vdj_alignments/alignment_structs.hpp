#pragma once

#include "gene_database.hpp"

// in our context, query is a reads, subject is a gene segment
struct AlignmentPositions {
    std::pair<size_t, size_t> query_pos;
    std::pair<size_t, size_t> subject_pos;

    AlignmentPositions(std::pair<size_t, size_t> new_query_pos,
              std::pair<size_t, size_t> new_subject_pos) :
            query_pos(new_query_pos),
            subject_pos(new_subject_pos) { }
};

std::ostream& operator<<(std::ostream &out, const AlignmentPositions& obj);

//-----------------------------------------------------------

struct IgGeneAlignmentPositions {
    AlignmentPositions alignment;
    const IgGene& ig_gene;

    IgGeneAlignmentPositions(AlignmentPositions new_alignment,
                   const IgGene& new_ig_gene) :
            alignment(new_alignment),
            ig_gene(new_ig_gene) { }
};

std::ostream& operator<<(std::ostream& out, const IgGeneAlignmentPositions& obj);

//-----------------------------------------------------------

//struct