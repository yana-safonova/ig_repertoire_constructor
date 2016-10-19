#include "poly_evolutionary_edge_constructor.hpp"

#include <annotation_utils/shm_comparator.hpp>

namespace antevolo {
    std::shared_ptr<BaseEvolutionaryEdge> PolySimpleEvolutionaryEdgeConstructor::ConstructEdge(
            const annotation_utils::AnnotatedClone &src_clone,
            const annotation_utils::AnnotatedClone &dst_clone,
            size_t src_num, size_t dst_num) const {

        if (annotation_utils::SHMComparator::SHMsAreEqual(src_clone.VSHMs(), dst_clone.VSHMs()))
            return std::make_shared( new UndirectedEvolutionaryEdge(src_clone, dst_clone, src_num, dst_num) );
        if (annotation_utils::SHMComparator::SHMs1AreNestedInSHMs2(src_clone.VSHMs(), dst_clone.VSHMs()))
            return std::make_shared( new DirectedEvolutionaryEdge(src_clone, dst_clone, src_num, dst_num) );
        if (annotation_utils::SHMComparator::GetNumberOfIntersections(src_clone.VSHMs(), dst_clone.VSHMs()) >=
            params_.min_num_intersected_v_shms)
            return std::make_shared( new IntersectedEvolutionaryEdge(src_clone, dst_clone, src_num, dst_num,
                                                                     params_.min_num_intersected_v_shms) );
        return std::make_shared( new BaseEvolutionaryEdge(src_clone, dst_clone, src_num, dst_num) );
    }

    std::shared_ptr<BaseEvolutionaryEdge> PolyVJEvolutionaryEdgeConstructor::ConstructEdge(
            const annotation_utils::AnnotatedClone &src_clone,
            const annotation_utils::AnnotatedClone &dst_clone,
            size_t src_num, size_t dst_num) const {

        if (annotation_utils::SHMComparator::SHMsAreEqual(src_clone.VSHMs(), dst_clone.VSHMs()) &&
            annotation_utils::SHMComparator::SHMsAreEqual(src_clone.JSHMs(), dst_clone.JSHMs()))
            return std::make_shared( new UndirectedEvolutionaryEdge(src_clone, dst_clone, src_num, dst_num) );
        if (annotation_utils::SHMComparator::SHMs1AreNestedInSHMs2(src_clone.VSHMs(), dst_clone.VSHMs()) &&
            annotation_utils::SHMComparator::SHMs1AreNestedInSHMs2(src_clone.JSHMs(), dst_clone.JSHMs()))
            return std::make_shared( new DirectedEvolutionaryEdge(src_clone, dst_clone, src_num, dst_num) );
        if (annotation_utils::SHMComparator::GetNumberOfIntersections(src_clone.VSHMs(), dst_clone.VSHMs()) >=
            params_.min_num_intersected_v_shms)
            return std::make_shared( new IntersectedEvolutionaryEdge(src_clone, dst_clone, src_num, dst_num,
                                                                     params_.min_num_intersected_v_shms) );
        return std::make_shared( new BaseEvolutionaryEdge(src_clone, dst_clone, src_num, dst_num) );
    }
}