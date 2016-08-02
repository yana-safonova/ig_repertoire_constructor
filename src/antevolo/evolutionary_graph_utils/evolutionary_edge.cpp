#include <annotation_utils/shm_comparator.hpp>
#include "evolutionary_edge.hpp"

namespace antevolo {
    //todo: refactor it
    template<typename T>
    size_t HammingDistance(T seq1, T seq2) {
        size_t dist = 0;
        size_t min_length = std::min<size_t>(seqan::length(seq1), seqan::length(seq2));
        for(size_t i = 0; i < min_length; i++) {
            if(seq1[i] != seq2[i])
                dist++;
        }
        return dist;
    }

    std::ostream& operator<<(std::ostream &out, const EvolutionaryEdgeType& edge_type) {
        if(edge_type == EvolutionaryEdgeType::DirectedEdgeType)
            out << "directed";
        else if(edge_type == EvolutionaryEdgeType::UndirectedEdgeType)
            out << "undirected";
        else if(edge_type == EvolutionaryEdgeType::IntersectedEdgeType)
            out << "intersected";
        else
            out << "unknown evolutionary edge";
        return out;
    }

    void EvolutionaryEdge::InitializeFields() {
        if(edge_type == EvolutionaryEdgeType::UnknownEdgeType)
            return;
        // todo: refactor it. now cdr3s can differ only by mismatches
        cdr3_distance = HammingDistance(src_clone->CDR3(), dst_clone->CDR3());
        if(edge_type == EvolutionaryEdgeType::UndirectedEdgeType) {
            num_added_v_shms = 0;
            num_intersected_v_shms = src_clone->VSHMs().size();
            num_added_j_shms = 0;
            num_intersected_j_shms = src_clone->JSHMs().size();
            num_added_shms = num_added_v_shms + num_added_j_shms;
            num_intersected_shms = num_intersected_v_shms + num_intersected_j_shms;
            weight = cdr3_distance;
        }
        else {
            if(edge_type == EvolutionaryEdgeType::DirectedEdgeType) {
                VERIFY_MSG(dst_clone->VSHMs().size() + dst_clone->JSHMs().size() >
                           src_clone->VSHMs().size() + src_clone->JSHMs().size(),
                           "# SHMs in destination clone (" <<
                                   dst_clone->VSHMs().size() + dst_clone->JSHMs().size() <<
                           ") does not exceed # SHMs in source clone (" <<
                                   src_clone->VSHMs().size() + src_clone->JSHMs().size() << ")");
                num_added_v_shms = dst_clone->VSHMs().size() - src_clone->VSHMs().size();
                num_intersected_v_shms = src_clone->VSHMs().size();
                num_added_j_shms = dst_clone->JSHMs().size() - src_clone->JSHMs().size();
                num_intersected_j_shms = src_clone->JSHMs().size();
                num_added_shms = num_added_v_shms + num_added_j_shms;
                num_intersected_shms = num_intersected_v_shms + num_intersected_j_shms;
                weight = cdr3_distance + num_added_shms;
            }
            else { // INTERSECTED CASE JSHMs have not been not added !!
                num_intersected_v_shms = annotation_utils::SHMComparator::GetNumberOfIntersections(src_clone->VSHMs(),
                                                                                                   dst_clone->VSHMs());
                num_added_v_shms = src_clone->VSHMs().size() + dst_clone->VSHMs().size() - 2 * num_intersected_v_shms;
                weight = cdr3_distance * intersected_edge_coeff_ + num_added_v_shms;
            }
        }
    }

    bool EvolutionaryEdge::IsSynonymous() const {
        return annotation_utils::SHMComparator::AddedSHMsAreSynonimous(src_clone->JSHMs(), dst_clone->JSHMs()) &&
               annotation_utils::SHMComparator::AddedSHMsAreSynonimous(src_clone->VSHMs(), dst_clone->VSHMs());
    }


    std::ostream& operator<<(std::ostream& out, const EvolutionaryEdge &edge) {
        out << "Type: " << edge.edge_type << ", # added V SHMs: " << edge.num_added_v_shms <<
                ", # shared SHMs: " << edge.num_intersected_v_shms <<
                ", CDR3 distance: " << edge.cdr3_distance << ", weight: " << edge.weight;
        return out;
    }
}