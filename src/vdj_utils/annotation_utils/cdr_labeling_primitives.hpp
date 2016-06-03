#pragma once

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

    struct CDRLabeling {
        CDRRange cdr1;
        CDRRange cdr2;
        CDRRange cdr3;

        CDRLabeling() : cdr1(), cdr2(), cdr3() { }

        CDRLabeling(CDRRange cdr1, CDRRange cdr2, CDRRange cdr3) :
                cdr1(cdr1), cdr2(cdr2), cdr3(cdr3) { }

        bool Empty() const { return cdr1.Empty() and cdr2.Empty() and cdr3.Empty(); }
    };
}