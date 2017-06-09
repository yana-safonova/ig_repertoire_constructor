#include "tree_annotated_shm.hpp"

namespace antevolo {
    bool operator<(const TreeSHM &left, const TreeSHM &right) {
        if(left.gene_pos != right.gene_pos)
            return left.gene_pos < right.gene_pos;
        if(left.src_triplet != right.src_triplet)
            return left.src_triplet < right.src_triplet;
        return left.dst_triplet < right.dst_triplet;
    }

    bool TreeSHM::operator==(const TreeSHM &obj) const {
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

    std::string GetTripletByCentralPos(seqan::Dna5String str, size_t pos, size_t orf) {
        size_t aa_pos = (pos - orf) / 3;
        std::stringstream ss;
        ss << str[aa_pos * 3 + orf] << str[aa_pos * 3 + orf + 1] << str[aa_pos * 3 + orf + 2];
        return ss.str();
    }

    TreeSHM TreeSHMCalculator::ComputeTreeSHM(annotation_utils::SHM shm, size_t src, size_t dst) const {
        TreeSHM tree_shm;
        // gene
        tree_shm.gene_pos = shm.gene_nucl_pos;
        tree_shm.gene_nucl = shm.gene_nucl;
        tree_shm.gene_aa = shm.gene_aa;
        tree_shm.gene_aa_pos = (shm.gene_nucl_pos - clone_set_[dst].VAlignment().subject().ORF()) / 3;
        //tree_shm.gene_name = clone_set_[dst].VAlignment().subject().name();
        // dst
        tree_shm.dst_pos = shm.read_nucl_pos;
        tree_shm.dst_nucl = shm.read_nucl;
        tree_shm.dst_aa = shm.read_aa;
        // src
        tree_shm.src_pos = clone_set_[src].VAlignment().QueryPositionBySubjectPosition(tree_shm.gene_pos);
        tree_shm.src_nucl = clone_set_[src].Read().seq[tree_shm.src_pos];
        tree_shm.src_aa = clone_set_[src].GetAminoAcidByNucleotidePos(tree_shm.src_pos);
        // triplets
        tree_shm.src_triplet = GetTripletByCentralPos(clone_set_[src].Read().seq, tree_shm.src_pos, clone_set_[src].ORF());
        tree_shm.dst_triplet = GetTripletByCentralPos(clone_set_[dst].Read().seq, tree_shm.dst_pos, clone_set_[dst].ORF());
        // region
        tree_shm.region = clone_set_[dst].GetRegionBySHM(shm);
        return tree_shm;
    }
}