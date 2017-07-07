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

        typedef typename std::vector<AnnotatedClone>::iterator iterator;
        typedef typename std::vector<AnnotatedClone>::const_iterator const_iterator;

        iterator       begin ()       { return annotated_clones_.begin (); }
        const_iterator begin () const { return annotated_clones_.begin (); }
        const_iterator cbegin() const { return annotated_clones_.cbegin(); }
        iterator       end   ()       { return annotated_clones_.end   (); }
        const_iterator end   () const { return annotated_clones_.end   (); }
        const_iterator cend  () const { return annotated_clones_.cend  (); }
        size_t size() const { return annotated_clones_.size(); }

        const AnnotatedClone& operator[](size_t index) const {
            VERIFY_MSG(index < size(), "Index " << index << " exceeds size of clone set (" << size() << ")");
            return annotated_clones_[index];
        }

        AnnotatedClone& operator[](size_t index) {
            VERIFY_MSG(index < size(), "Index " << index << " exceeds size of clone set (" << size() << ")");
            return annotated_clones_[index];
        }
    };

    typedef AnnotatedCloneSet<AnnotatedClone> CDRAnnotatedCloneSet;
}