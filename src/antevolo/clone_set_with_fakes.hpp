#pragma once

#include "../vdj_utils/annotation_utils/annotated_clone.hpp"
#include "../vdj_utils/annotation_utils/annotated_clone_set.hpp"

namespace antevolo {
    class CloneSetWithFakes : public annotation_utils::AnnotatedCloneSet<annotation_utils::AnnotatedClone> {
        size_t first_fake_;

    public:
        explicit CloneSetWithFakes(const annotation_utils::AnnotatedCloneSet<annotation_utils::AnnotatedClone>& clone_set) :
                first_fake_(0) {
            for (auto it = clone_set.cbegin(); it != clone_set.cend(); it++) {
                AddClone(*it);
                ++first_fake_;
            }
        }
        bool IsFake(size_t num) {
            return num >= first_fake_;
        }

        void PopBack() {
            annotated_clones_[]
        }
    };
}
