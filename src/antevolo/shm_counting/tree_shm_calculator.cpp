#include "tree_shm_calculator.hpp"

namespace antevolo {
    TreeSHM TreeSHMCalculator::ComputeTreeSHMByUsualSHM(annotation_utils::SHM shm, size_t src, size_t dst) const {
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

    void UniqueSHMCalculator::AddAddedSHMs(std::vector <annotation_utils::SHM> &shms, size_t src, size_t dst) {
        for(auto it = shms.begin(); it != shms.end(); it++) {
            shm_map_.AddSHM(tree_shm_calc_.ComputeTreeSHMByUsualSHM(*it, src, dst), src, dst);
        }
    }

    void UniqueSHMCalculator::AddNestedSHMs(size_t src, size_t dst) {
        auto v_shms = annotation_utils::SHMComparator::GetAddedSHMs(clone_set_[src].VSHMs(), clone_set_[dst].VSHMs());
        AddAddedSHMs(v_shms, src, dst);
        auto j_shms = annotation_utils::SHMComparator::GetAddedSHMs(clone_set_[src].JSHMs(), clone_set_[dst].JSHMs());
        AddAddedSHMs(j_shms, src, dst);
    }

    TreeSHM UniqueSHMCalculator::ComputeSHMInCDR3(size_t src, size_t dst, size_t cdr3_pos) const {
        size_t src_pos = clone_set_[src].CDR3Range().start_pos + cdr3_pos;
        size_t dst_pos = clone_set_[dst].CDR3Range().start_pos + cdr3_pos;
        TreeSHM tree_shm;
        tree_shm.gene_pos = cdr3_pos;
        tree_shm.src_pos = src_pos;
        tree_shm.dst_pos = dst_pos;
        tree_shm.gene_nucl = '-';
        tree_shm.src_nucl = clone_set_[src].Read().seq[src_pos];
        tree_shm.dst_nucl = clone_set_[dst].Read().seq[dst_pos];
        tree_shm.gene_aa = '-';
        tree_shm.src_aa = clone_set_[src].GetAminoAcidByNucleotidePos(src_pos);
        tree_shm.dst_aa = clone_set_[dst].GetAminoAcidByNucleotidePos(dst_pos);
        tree_shm.src_triplet = GetTripletByCentralPos(clone_set_[src].Read().seq, src_pos, clone_set_[src].ORF());
        tree_shm.dst_triplet = GetTripletByCentralPos(clone_set_[dst].Read().seq, dst_pos, clone_set_[dst].ORF());
        tree_shm.region = annotation_utils::StructuralRegion::CDR3;
        return tree_shm;
    }

    // todo: add in config
    void UniqueSHMCalculator::AddCDR3SHMs(size_t src, size_t dst) {
        auto cdr3_src = clone_set_[src].CDR3();
        auto cdr3_dst = clone_set_[dst].CDR3();
        if(seqan::length(cdr3_dst) != seqan::length(cdr3_src)) {
            WARN("CDR3s of clones #" << src << "(" << cdr3_src << ") & " << dst << "(" << cdr3_dst <<
                                     ") have different lengths");
            return;
        }
        for(size_t i = 0; i < seqan::length(cdr3_dst); i++) {
            if(cdr3_src[i] != cdr3_dst[i]) {
                shm_map_.AddSHM(ComputeSHMInCDR3(src, dst, i), src, dst);
            }
        }
    }

    // todo: add as separate filter
    bool UniqueSHMCalculator::FilterEdge(size_t src, size_t dst) const {
        return src >= clone_set_.size() or dst >= clone_set_.size();
    }

    //TreeSHMMap UniqueSHMCalculator::ComputeSHMMap() {
    //    for(auto it = tree_.cbegin(); it != tree_.cend(); it++) {
    //        auto src = (*it)->SrcNum();
    //        auto dst = (*it)->DstNum();
    //        if(FilterEdge(src, dst))
    //            continue;
    //        if((*it)->IsDirected())
    //            AddNestedSHMs(src, dst);
    //        AddCDR3SHMs(src, dst);
    //    }
    //    return shm_map_;
    //}

    void UniqueSHMCalculator::AddSHMsFromEdge(const EvolutionaryEdgePtr edge_ptr) {
        AddSHMsFromEdge(edge_ptr->SrcNum(), edge_ptr->DstNum());
    }

    void UniqueSHMCalculator::AddSHMsFromEdge(size_t src_id, size_t dst_id) {
        //if(edge_ptr->IsDirected())
        AddNestedSHMs(src_id, dst_id);
        AddNestedSHMs(dst_id, src_id);
        AddCDR3SHMs(src_id, dst_id);
    }
}