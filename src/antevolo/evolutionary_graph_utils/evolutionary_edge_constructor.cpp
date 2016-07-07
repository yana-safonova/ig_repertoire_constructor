#include "evolutionary_edge_constructor.hpp"

#include <annotation_utils/shm_comparator.hpp>

namespace antevolo {
    EvolutionaryEdge SimpleEvolutionaryEdgeConstructor::ConstructEdge(const annotation_utils::AnnotatedClone &src_clone,
                                                                      const annotation_utils::AnnotatedClone &dst_clone) const {
        if(annotation_utils::SHMComparator::SHMsAreEqual(src_clone.VSHMs(), dst_clone.VSHMs()))
            return EvolutionaryEdge(EvolutionaryEdgeType::UndirectedEdgeType, src_clone, dst_clone,
                                    params_.intersected_edge_coeff);
        if(annotation_utils::SHMComparator::SHMs1AreNestedInSHMs2(src_clone.VSHMs(), dst_clone.VSHMs()))
            return EvolutionaryEdge(EvolutionaryEdgeType::DirectedEdgeType, src_clone, dst_clone,
                                    params_.intersected_edge_coeff);
        if(annotation_utils::SHMComparator::GetNumberOfIntersections(src_clone.VSHMs(), dst_clone.VSHMs()) >=
                params_.min_num_intersected_v_shms)
            return EvolutionaryEdge(EvolutionaryEdgeType::IntersectedEdgeType, src_clone, dst_clone,
                                    params_.intersected_edge_coeff);
        return EvolutionaryEdge(EvolutionaryEdgeType::UnknownEdgeType, src_clone, dst_clone,
                                params_.intersected_edge_coeff);
    }
}