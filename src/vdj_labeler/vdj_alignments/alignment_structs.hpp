#pragma once

#include <seqan/align.h>

#include "gene_database.hpp"
#include "../fastq_read_archive.hpp"

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
    IgGenePtr ig_gene;
    ReadPtr read;

    IgGeneAlignmentPositions(AlignmentPositions new_alignment,
                             IgGenePtr new_ig_gene,
                             ReadPtr new_read) :
            alignment(new_alignment),
            ig_gene(new_ig_gene),
            read(new_read) { }
};

std::ostream& operator<<(std::ostream& out, const IgGeneAlignmentPositions& obj);

//-----------------------------------------------------------

struct IgGeneAlignment {
    IgGeneAlignmentPositions positions;
    seqan::Align<Dna5String> alignment;

    IgGeneAlignment(IgGeneAlignmentPositions new_positions,
                    seqan::Align<Dna5String> new_alignment) :
            positions(new_positions),
            alignment(new_alignment) { }
};

std::ostream& operator<<(std::ostream &out, const IgGeneAlignment& ig_gene_alignment);

typedef std::shared_ptr<IgGeneAlignment> IgGeneAlignmentPtr;