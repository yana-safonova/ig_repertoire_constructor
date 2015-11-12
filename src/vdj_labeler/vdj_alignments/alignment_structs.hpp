#pragma once

#include "gene_database.hpp"

struct Alignment {
    std::pair<size_t, size_t> query_pos;
    std::pair<size_t, size_t> subject_pos;

    Alignment(std::pair<size_t, size_t> new_query_pos,
              std::pair<size_t, size_t> new_subject_pos) :
            query_pos(new_query_pos),
            subject_pos(new_subject_pos) { }
};

struct IgGeneAlignment {
    Alignment alignment;
    const IgGene& ig_gene;

    IgGeneAlignment(Alignment new_alignment,
                   const IgGene& new_ig_gene) :
            alignment(new_alignment),
            ig_gene(new_ig_gene) { }
};