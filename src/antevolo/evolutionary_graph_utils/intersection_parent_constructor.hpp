#pragma once

#include <annotation_utils/annotated_clone.hpp>

namespace antevolo {

    class IntersectionParentConstructor {
    public:
        static annotation_utils::AnnotatedClone ReconstructParent(
                const std::shared_ptr<annotation_utils::AnnotatedClone>& clone1,
                const std::shared_ptr<annotation_utils::AnnotatedClone>& clone2);
    };

}
