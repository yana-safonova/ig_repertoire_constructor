#include "alignment_structs.hpp"

using namespace std;

ostream& operator<<(ostream &out, const AlignmentPositions& obj) {
    out << "Query pos: " << obj.query_pos.first << " - " << obj.query_pos.second << ". Subject pos: " <<
            obj.subject_pos.first << " - " << obj.subject_pos.second;
    return out;
}

ostream& operator<<(ostream& out, const IgGeneAlignmentPositions& obj) {
    out << obj.ig_gene << endl;
    out << obj.alignment;
    return out;
}