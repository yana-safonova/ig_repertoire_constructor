#include "evolutionary_edge_constructor.hpp"

#include <annotation_utils/shm_comparator.hpp>

namespace antevolo {
//    std::shared_ptr<BaseEvolutionaryEdge> SimpleEvolutionaryEdgeConstructor::ConstructEdge(
//            const annotation_utils::AnnotatedClone &src_clone,
//            const annotation_utils::AnnotatedClone &dst_clone,
//            size_t src_num, size_t dst_num) const {
//
//        if (annotation_utils::SHMComparator::SHMsAreEqual(src_clone.VSHMs(), dst_clone.VSHMs()))
//            return std::shared_ptr<BaseEvolutionaryEdge>( new UndirectedEvolutionaryEdge(src_clone, dst_clone,
//                                                                                         src_num, dst_num) );
//        if (annotation_utils::SHMComparator::SHMs1AreNestedInSHMs2(src_clone.VSHMs(), dst_clone.VSHMs()))
//            return std::shared_ptr<BaseEvolutionaryEdge>( new DirectedEvolutionaryEdge(src_clone, dst_clone,
//                                                                                       src_num, dst_num) );
//        if (annotation_utils::SHMComparator::GetNumberOfIntersections(src_clone.VSHMs(), dst_clone.VSHMs()) >=
//            params_.min_num_intersected_v_shms)
//            return std::shared_ptr<BaseEvolutionaryEdge>( new IntersectedEvolutionaryEdge(src_clone, dst_clone,
//                                                                                          src_num, dst_num) );
//        return std::shared_ptr<BaseEvolutionaryEdge>( new BaseEvolutionaryEdge(src_clone, dst_clone,
//                                                                               src_num, dst_num) );
//    }

    std::shared_ptr<BaseEvolutionaryEdge> VJEvolutionaryEdgeConstructor::ConstructEdge(
            const annotation_utils::AnnotatedClone &src_clone,
            const annotation_utils::AnnotatedClone &dst_clone,
            size_t src_num, size_t dst_num) const {
        if (!(annotation_utils::SHMComparator::SHMsInsertionBlocksAreEqual(src_clone.VSHMs(),
                                                                           dst_clone.VSHMs()) &&
              annotation_utils::SHMComparator::SHMsInsertionBlocksAreEqual(src_clone.JSHMs(),
                                                                           dst_clone.JSHMs()))   ||
             src_clone.CDR3Range().length() != dst_clone.CDR3Range().length())
            {
            return std::shared_ptr<BaseEvolutionaryEdge>( new BaseEvolutionaryEdge(src_clone, dst_clone,
                                                                                   src_num, dst_num) );
        }
        // undirected
        if (annotation_utils::SHMComparator::SHMsAreEqual(src_clone.VSHMs(), dst_clone.VSHMs()) &&
            annotation_utils::SHMComparator::SHMsAreEqual(src_clone.JSHMs(), dst_clone.JSHMs()))
            return std::shared_ptr<BaseEvolutionaryEdge>( new UndirectedEvolutionaryEdge(src_clone, dst_clone,
                                                                                         src_num, dst_num) );
        // directed
        if (annotation_utils::SHMComparator::SHMs1AreNestedInSHMs2(src_clone.VSHMs(), dst_clone.VSHMs()) &&
            annotation_utils::SHMComparator::SHMs1AreNestedInSHMs2(src_clone.JSHMs(), dst_clone.JSHMs()))
            return std::shared_ptr<BaseEvolutionaryEdge>( new DirectedEvolutionaryEdge(src_clone, dst_clone,
                                                                                       src_num, dst_num) );
        // reverse directed
        if (annotation_utils::SHMComparator::SHMs1AreNestedInSHMs2(dst_clone.VSHMs(), src_clone.VSHMs()) &&
            annotation_utils::SHMComparator::SHMs1AreNestedInSHMs2(dst_clone.JSHMs(), src_clone.JSHMs()))
            return std::shared_ptr<BaseEvolutionaryEdge>( new ReverseDirectedEvolutionaryEdge(src_clone,
                                                                                              dst_clone,
                                                                                              src_num,
                                                                                              dst_num) );
        // intersected
        // V
        size_t num_intersected_v_shms = annotation_utils::SHMComparator::GetNumberOfIntersections(src_clone.VSHMs(),
                                                                                                  dst_clone.VSHMs());
        size_t num_individual_v_shms = src_clone.VSHMs().size() + dst_clone.VSHMs().size() - 2 * num_intersected_v_shms;

        // J
        size_t num_intersected_j_shms = annotation_utils::SHMComparator::GetNumberOfIntersections(src_clone.JSHMs(),
                                                                                           dst_clone.JSHMs());
        size_t num_individual_j_shms = src_clone.JSHMs().size() + dst_clone.JSHMs().size() - 2 * num_intersected_j_shms;
        // sum
        size_t num_individual_shms = num_individual_v_shms + num_individual_v_shms;
        size_t num_shared_shms = num_intersected_v_shms + num_intersected_j_shms;
        if (num_shared_shms >= params_.min_num_intersected_v_shms &&
            num_individual_shms <= num_shared_shms * params_.min_num_intersected_v_shms)
            return std::shared_ptr<BaseEvolutionaryEdge>( new IntersectedEvolutionaryEdge(src_clone, dst_clone,
                                                                                          src_num, dst_num,
                                                                                          num_individual_v_shms,
                                                                                          num_intersected_v_shms,
                                                                                          num_individual_j_shms,
                                                                                          num_intersected_j_shms) );
        return std::shared_ptr<BaseEvolutionaryEdge>( new BaseEvolutionaryEdge(src_clone, dst_clone,
                                                                               src_num, dst_num) );
    }
}
