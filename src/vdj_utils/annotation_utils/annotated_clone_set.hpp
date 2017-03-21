#pragma once

#include <vector>
#include <verify.hpp>

#include "annotated_clone.hpp"

namespace annotation_utils {
    template<typename AnnotatedClone>
    class AnnotatedCloneSet {
        std::vector<AnnotatedClone> annotated_clones_;

    public:
        void AddClone(AnnotatedClone clone) {
            annotated_clones_.push_back(clone);
        }

        typedef typename std::vector<AnnotatedClone>::const_iterator AnnotatedCloneIterator;

        AnnotatedCloneIterator cbegin() const { return annotated_clones_.cbegin(); }

        AnnotatedCloneIterator cend() const { return annotated_clones_.cend(); }

        size_t size() const { return annotated_clones_.size(); }

        const AnnotatedClone& operator[](size_t index) const {
            VERIFY_MSG(index < size(), "Index " << index << " exceeds size of clone set");
            return annotated_clones_[index];
        }
    };

    typedef AnnotatedCloneSet<AnnotatedClone> CDRAnnotatedCloneSet;
}