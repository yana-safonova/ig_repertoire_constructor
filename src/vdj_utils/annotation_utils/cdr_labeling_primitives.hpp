#pragma once

#include <ostream>

namespace annotation_utils {
    // start and end positions are inclusive
    struct CDRRange {
        size_t start_pos;
        size_t end_pos;

        CDRRange() :
                start_pos(size_t(-1)),
                end_pos(size_t(-1)) { }

        CDRRange(size_t start_pos, size_t end_pos) :
                start_pos(start_pos),
                end_pos(end_pos) { }

        bool Valid() const {
            if(start_pos == size_t(-1) and end_pos == size_t(-1))
                return false;
            if(start_pos != size_t(-1) and end_pos != size_t(-1))
                return start_pos < end_pos;
            return true;
        }

        bool Empty() const { return !Valid(); }

        bool Full() const { return start_pos != size_t(-1) and end_pos != size_t(-1); }

        size_t length() const;
    };

    std::ostream& operator<<(std::ostream& out, const CDRRange &obj);

    struct CDRLabeling {
        CDRRange cdr1;
        CDRRange cdr2;
        CDRRange cdr3;
        unsigned orf;

        CDRLabeling() : cdr1(), cdr2(), cdr3(), orf(0) { }

        CDRLabeling(CDRRange cdr1, CDRRange cdr2, CDRRange cdr3, unsigned orf = 0) :
                cdr1(cdr1), cdr2(cdr2), cdr3(cdr3), orf(orf) { }

        bool Empty() const { return cdr1.Empty() and cdr2.Empty() and cdr3.Empty(); }
    };
}