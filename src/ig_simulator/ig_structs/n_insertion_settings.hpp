#pragma once

#include "../include_me.hpp"

class HC_NInsertionSettings {
    string vd_insertion_;
    string dj_insertion_;

public:
    HC_NInsertionSettings() :
            vd_insertion_(),
            dj_insertion_() { }

    HC_NInsertionSettings(string vd_insertion, string dj_insertion) :
            vd_insertion_(vd_insertion),
            dj_insertion_(dj_insertion) { }

    string VD_Insertion() const { return vd_insertion_; }

    string DJ_Insertion() const { return dj_insertion_; }
};

ostream& operator<<(ostream &out, const HC_NInsertionSettings &obj) {
    out << "VD insertion: " << obj.VD_Insertion() << ", len: " << obj.VD_Insertion().size() << endl;
    out << "DJ insertion: " << obj.DJ_Insertion() << ", len: " << obj.DJ_Insertion().size() << endl;
    return out;
}

// ----------------------------------------------------------------------------

class LC_NInsertionSettings {
    string vj_insertion_;

public:
    LC_NInsertionSettings() :
            vj_insertion_() { }

    LC_NInsertionSettings(string vj_insertion) :
            vj_insertion_(vj_insertion) { }

    string VJ_Insertion() const { return vj_insertion_; }
};

ostream& operator<<(ostream &out, const LC_NInsertionSettings &obj) {
    out << "VJ insertion: " << obj.VJ_Insertion() << ", len: " << obj.VJ_Insertion().size() << endl;
    return out;
}