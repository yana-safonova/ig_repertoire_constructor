#pragma once

#include "../include_me.hpp"

class CDRLabel {
    size_t start_;
    size_t end_;

public:
    CDRLabel() :
        start_(0), end_(0) { }

    CDRLabel(size_t start, size_t end) :
            start_(start), end_(end) { }

    size_t Start() const { return start_; }

    size_t End() const { return end_; }
};

ostream& operator<<(ostream &out, const CDRLabel &label) {
    out << label.Start() << " - " << label.End();
    return out;
}

// ----------------------------------------------------------

class CDRSettings {
    CDRLabel cdr1_;
    CDRLabel cdr2_;
    CDRLabel cdr3_;

public:
    CDRSettings() :
        cdr1_(), cdr2_(), cdr3_() { }

    CDRSettings(CDRLabel cdr1, CDRLabel cdr2, CDRLabel cdr3):
            cdr1_(cdr1),
            cdr2_(cdr2),
            cdr3_(cdr3) { }

    CDRLabel CDR1() const { return cdr1_; }

    CDRLabel CDR2() const { return cdr2_; }

    CDRLabel CDR3() const { return cdr3_; }
};

ostream& operator<<(ostream &out, const CDRSettings &cdr_settings) {
    out << "CDR1: " << cdr_settings.CDR1() << endl;
    out << "CDR2: " << cdr_settings.CDR2() << endl;
    out << "CDR3: " << cdr_settings.CDR3() << endl;
    return out;
}