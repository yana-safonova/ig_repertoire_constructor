#pragma once

#include "annotation_utils/annotated_clone_set.hpp"

namespace antevolo {
    struct TreeSHM {
        size_t gene_pos;
        size_t src_pos;
        size_t dst_pos;

        char gene_nucl;
        char src_nucl;
        char dst_nucl;

        char gene_aa;
        char src_aa;
        char dst_aa;

        size_t gene_aa_pos;

        std::string src_triplet;
        std::string dst_triplet;

        annotation_utils::StructuralRegion region;

    public:
        TreeSHM() { }

        TreeSHM(size_t gene_pos, size_t src_pos, size_t dst_pos,
                char gene_nucl, char src_nucl, char dst_nucl,
                char gene_aa, char src_aa, char dst_aa, size_t gene_aa_pos,
                std::string src_triplet, std::string dst_triplet,
                annotation_utils::StructuralRegion region) : gene_pos(gene_pos), src_pos(src_pos), dst_pos(dst_pos),
                                                             gene_nucl(gene_nucl), src_nucl(src_nucl), dst_nucl(dst_nucl),
                                                             gene_aa(gene_aa), src_aa(src_aa), dst_aa(dst_aa),
                                                             gene_aa_pos(gene_aa_pos),
                                                             src_triplet(src_triplet), dst_triplet(dst_triplet),
                                                             region(region) { }

        bool Synonymous() const { return src_aa == dst_aa; }

        bool ToStopCodon() const { return dst_aa == '*'; }

        bool operator==(const TreeSHM &obj) const;
    };

    bool operator<(const TreeSHM &left, const TreeSHM &right);

    std::ostream& operator<<(std::ostream& out, const TreeSHM &tree_shm);

    std::string GetTripletByCentralPos(seqan::Dna5String str, size_t pos, size_t orf);
}