#pragma once

#include "../vdj_utils/annotation_utils/annotated_clone.hpp"
#include "../vdj_utils/annotation_utils/annotated_clone_set.hpp"

namespace antevolo {
    class CloneSetWithFakes {
        const annotation_utils::CDRAnnotatedCloneSet& original_clone_set_;
        annotation_utils::CDRAnnotatedCloneSet fakes_clone_set_;
        size_t first_fake_;

    public:
        explicit CloneSetWithFakes(const annotation_utils::CDRAnnotatedCloneSet& original_clone_set) :
                original_clone_set_(original_clone_set),
                first_fake_(original_clone_set_.size()) {}

        bool IsFake(size_t i) const {
            return i >= first_fake_;
        }
        const annotation_utils::AnnotatedClone& operator[](size_t i) const {
            if (IsFake(i)) {
                return fakes_clone_set_[i - first_fake_];
            }
            return original_clone_set_[i];
        }
        annotation_utils::AnnotatedClone& operator[](size_t i) {
            if (IsFake(i)) {
                return fakes_clone_set_[i - first_fake_];
            }
            return const_cast<annotation_utils::AnnotatedClone&>(original_clone_set_[i]);
        }


        void AddClone(const annotation_utils::AnnotatedClone& clone) {
            fakes_clone_set_.AddClone(clone);
        }
        size_t size() const { return first_fake_+ fakes_clone_set_.size(); }
        const annotation_utils::CDRAnnotatedCloneSet& GetOriginalCloneSet() const { return original_clone_set_; }
    };

    typedef std::shared_ptr<CloneSetWithFakes>  CloneSetWithFakesPtr;
}
