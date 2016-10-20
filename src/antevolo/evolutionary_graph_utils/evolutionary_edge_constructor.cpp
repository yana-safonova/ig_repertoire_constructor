#include "evolutionary_edge_constructor.hpp"

#include <annotation_utils/shm_comparator.hpp>

namespace antevolo {
    std::shared_ptr<BaseEvolutionaryEdge> SimpleEvolutionaryEdgeConstructor::ConstructEdge(
            const annotation_utils::AnnotatedClone &src_clone,
            const annotation_utils::AnnotatedClone &dst_clone,
            size_t src_num, size_t dst_num) const {
        
        if (annotation_utils::SHMComparator::SHMsAreEqual(src_clone.VSHMs(), dst_clone.VSHMs()))
            return std::shared_ptr<BaseEvolutionaryEdge>( new UndirectedEvolutionaryEdge(src_clone, dst_clone, src_num, dst_num) );
        if (annotation_utils::SHMComparator::SHMs1AreNestedInSHMs2(src_clone.VSHMs(), dst_clone.VSHMs()))
            return std::shared_ptr<BaseEvolutionaryEdge>( new DirectedEvolutionaryEdge(src_clone, dst_clone, src_num, dst_num) );
        if (annotation_utils::SHMComparator::GetNumberOfIntersections(src_clone.VSHMs(), dst_clone.VSHMs()) >=
            params_.min_num_intersected_v_shms)
            return std::shared_ptr<BaseEvolutionaryEdge>( new IntersectedEvolutionaryEdge(src_clone, dst_clone, src_num, dst_num) );
        return std::shared_ptr<BaseEvolutionaryEdge>( new BaseEvolutionaryEdge(src_clone, dst_clone, src_num, dst_num) );
    }

    std::shared_ptr<BaseEvolutionaryEdge> VJEvolutionaryEdgeConstructor::ConstructEdge(
            const annotation_utils::AnnotatedClone &src_clone,
            const annotation_utils::AnnotatedClone &dst_clone,
            size_t src_num, size_t dst_num) const {

        if (annotation_utils::SHMComparator::SHMsAreEqual(src_clone.VSHMs(), dst_clone.VSHMs()) &&
            annotation_utils::SHMComparator::SHMsAreEqual(src_clone.JSHMs(), dst_clone.JSHMs()))
            return std::shared_ptr<BaseEvolutionaryEdge>( new UndirectedEvolutionaryEdge(src_clone, dst_clone, src_num, dst_num) );
        if (annotation_utils::SHMComparator::SHMs1AreNestedInSHMs2(src_clone.VSHMs(), dst_clone.VSHMs()) &&
            annotation_utils::SHMComparator::SHMs1AreNestedInSHMs2(src_clone.JSHMs(), dst_clone.JSHMs()))
            return std::shared_ptr<BaseEvolutionaryEdge>( new DirectedEvolutionaryEdge(src_clone, dst_clone, src_num, dst_num) );
        if (annotation_utils::SHMComparator::GetNumberOfIntersections(src_clone.VSHMs(), dst_clone.VSHMs()) >=
            params_.min_num_intersected_v_shms)
            return std::shared_ptr<BaseEvolutionaryEdge>( new IntersectedEvolutionaryEdge(src_clone, dst_clone, src_num, dst_num) );
        return std::shared_ptr<BaseEvolutionaryEdge>( new BaseEvolutionaryEdge(src_clone, dst_clone, src_num, dst_num) );
    }
}