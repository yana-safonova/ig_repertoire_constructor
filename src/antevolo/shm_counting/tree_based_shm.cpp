#include "tree_based_shm.hpp"

namespace antevolo {
    std::string GetTripletByCentralPos(seqan::Dna5String str, size_t pos, size_t orf) {
        size_t aa_pos = (pos - orf) / 3;
        std::stringstream ss;
        ss << str[aa_pos * 3 + orf] << str[aa_pos * 3 + orf + 1] << str[aa_pos * 3 + orf + 2];
        return ss.str();
    }

    bool operator<(const TreeSHM &left, const TreeSHM &right) {
        if(left.region != right.region)
            return left.region < right.region;
        if(left.gene_pos != right.gene_pos)
            return left.gene_pos < right.gene_pos;
        if(left.src_triplet != right.src_triplet)
            return left.src_triplet < right.src_triplet;
        return left.dst_triplet < right.dst_triplet;
    }

    bool TreeSHM::operator==(const TreeSHM &obj) const {
        if(region != obj.region)
            return false;
        if(gene_pos != obj.gene_pos)
            return false;
        if(src_triplet != obj.src_triplet)
            return false;
        return dst_triplet == obj.dst_triplet;
    }

    std::ostream& operator<<(std::ostream& out, const TreeSHM &tree_shm) {
        out << tree_shm.gene_pos << "(" << tree_shm.src_pos << "->" << tree_shm.dst_pos << "): " <<
            tree_shm.src_triplet << "->" << tree_shm.dst_triplet << ", " << tree_shm.src_aa << "->" <<
            tree_shm.dst_aa;
        return out;
    }
}