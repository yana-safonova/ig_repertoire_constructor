#include "alignment_structs.hpp"

using namespace std;

ostream& operator<<(ostream &out, const Alignment& obj) {
    out << "Query pos: " << obj.query_pos.first << " - " << obj.query_pos.second << ". Subject pos: " <<
            obj.subject_pos.first << " - " << obj.subject_pos.second;
    return out;
}

ostream& operator<<(ostream& out, const IgGeneAlignment& obj) {
    out << obj.ig_gene << endl;
    out << obj.alignment;
    return out;
}