#pragma once

#include "../vdj_utils/annotation_utils/annotated_clone.hpp"
#include "../vdj_utils/annotation_utils/annotated_clone_set.hpp"

namespace antevolo {
    class CloneSetWithFakes {
        const annotation_utils::CDRAnnotatedCloneSet& original_clone_set_;
        size_t first_fake_;
        annotation_utils::CDRAnnotatedCloneSet fakes_clone_set_;
    public:
        explicit CloneSetWithFakes(const annotation_utils::CDRAnnotatedCloneSet& original_clone_set) :
                original_clone_set_(original_clone_set),
                first_fake_(original_clone_set_.size()) {}

        bool IsFake(size_t num) const {
            return num >= first_fake_;
        }
        size_t GetFirstFakeIndex() const {
            return first_fake_;
        }
        const annotation_utils::AnnotatedClone& operator[](size_t i) const {
            if (i < first_fake_) {
                return original_clone_set_[i];
            }
            return fakes_clone_set_[i - first_fake_];
        }
        size_t size() const {
            return first_fake_ + fakes_clone_set_.size();
        }
        void AddClone(const annotation_utils::AnnotatedClone& clone) {
            fakes_clone_set_.AddClone(clone);
        }
    };
}
