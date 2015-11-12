#include "vj_alignment_info.hpp"
#include "

void VJAlignmentInfo::ParseLine(std::string line) {

}

VJAlignmentInfo::VJAlignmentInfo(std::string filename) {
    ifstream input(filename);
    if(!input.good()) { }
    while(!input.eof()) {
        std::string tmp;
        getline(input, tmp);
    }
}