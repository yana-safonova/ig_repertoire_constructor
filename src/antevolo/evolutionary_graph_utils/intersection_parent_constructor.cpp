#include "intersection_parent_constructor.hpp"

namespace antevolo {
    static annotation_utils::AnnotatedClone IntersectionParentConstructor::ReconstructParent(
            const std::shared_ptr<annotation_utils::AnnotatedClone>& clone1,
            const std::shared_ptr<annotation_utils::AnnotatedClone>& clone2) {
        //todo: implement. now it just returns clone1!!
        return annotation_utils::AnnotatedClone(clone1);
    }
}