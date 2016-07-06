#include "alignment_positions.hpp"

namespace alignment_utils {

std::ostream& operator<<(std::ostream &out, const AlignmentPositions &obj) {
    out << "Query pos: "     << obj.query_pos.first   << " - " << obj.query_pos.second
        << ". Subject pos: " << obj.subject_pos.first << " - " << obj.subject_pos.second;
    return out;
}

std::ostream& operator<<(std::ostream& out, const ImmuneGeneAlignmentPositions& obj) {
    out << obj.AlignmentPositions();
    return out;
}

}