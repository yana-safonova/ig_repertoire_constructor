#pragma once

#include "../include_me.hpp"

class HC_PInsertionSettings {
    string v_end_;
    string d_start_;
    string d_end_;
    string j_start_;

public:
    HC_PInsertionSettings() :
            v_end_(),
            d_start_(),
            d_end_(),
            j_start_() {
    }

    HC_PInsertionSettings(string v_end, string d_start, string d_end, string j_start) :
            v_end_(v_end),
            d_start_(d_start),
            d_end_(d_end),
            j_start_(j_start) { }

    string VEnd() const {
        return v_end_;
    }

    string DStart() const {
        return d_start_;
    }

    string DEnd() const {
        return d_end_;
    }

    string JStart() const {
        return j_start_;
    }
};

ostream& operator<<(ostream &out, const HC_PInsertionSettings &obj) {
    out << "V end (palindrome nt): " << obj.VEnd() << endl;
    out << "D start (palindrome nt): " << obj.DStart() << endl;
    out << "D end (palindrome nt): " << obj.DEnd() << endl;
    out << "J start (palindrome nt): " << obj.JStart() << endl;
    return out;
}

// ----------------------------------------------------------------------------

class LC_PInsertionSettings {
    string v_end_;
    string j_start_;

public:
    LC_PInsertionSettings() :
            v_end_(),
            j_start_() {
    }

    LC_PInsertionSettings(string v_end, string j_start) :
            v_end_(v_end),
            j_start_(j_start) { }

    string VEnd() const {
        return v_end_;
    }

    string JStart() const {
        return j_start_;
    }
};

ostream& operator<<(ostream &out, const LC_PInsertionSettings &obj) {
    out << "V end (palindrome nt): " << obj.VEnd() << endl;
    out << "J start (palindrome nt): " << obj.JStart() << endl;
    return out;
}