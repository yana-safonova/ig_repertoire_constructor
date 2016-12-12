#pragma once

#include "../include_me.hpp"

class HC_RemovingSettings {
    size_t v_end_len_;
    size_t d_start_len_;
    size_t d_end_len_;
    size_t j_start_len_;

public:
    HC_RemovingSettings():
            v_end_len_(),
            d_start_len_(),
            d_end_len_(),
            j_start_len_() { }

    HC_RemovingSettings(size_t v_end_len, size_t d_start_len, size_t d_end_len, size_t j_start_len) :
            v_end_len_(v_end_len),
            d_start_len_(d_start_len),
            d_end_len_(d_end_len),
            j_start_len_(j_start_len) { }

    size_t VEndLen() const { return v_end_len_; }
    size_t DStartLen() const { return d_start_len_; }
    size_t DEndLen() const { return d_end_len_; }
    size_t JStartLen() const { return j_start_len_; }
};

ostream& operator<<(ostream &out, const HC_RemovingSettings &obj) {
    out << "V end (# removed): " << obj.VEndLen() << endl;
    out << "D start (# removed): " << obj.DStartLen() << endl;
    out << "D end (# removed): " << obj.DEndLen() << endl;
    out << "J start (# removed): " << obj.JStartLen() << endl;
    return out;
}

// ----------------------------------------------------------------------------

class LC_RemovingSettings {
    size_t v_end_len_;
    size_t j_start_len_;

public:
    LC_RemovingSettings() :
            v_end_len_(),
            j_start_len_() { }

    LC_RemovingSettings(size_t v_end_len, size_t j_start_len) :
            v_end_len_(v_end_len),
            j_start_len_(j_start_len) { }

    size_t VEndLen() const { return v_end_len_; }
    size_t JStartLen() const { return j_start_len_; }
};

ostream& operator<<(ostream &out, const LC_RemovingSettings &obj) {
    out << "V end (# removed): " << obj.VEndLen() << endl;
    out << "J start (# removed): " << obj.JStartLen() << endl;
    return out;
}