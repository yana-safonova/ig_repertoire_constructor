#pragma once

#include "gene_segment_alignment.hpp"

class VJAlignmentInfo {
    std::vector<IgGeneAlignment> v_segments_;
    std::vector<IgGeneAlignment> j_segments_;
    std::vector<std::string> read_names_;

    void ParseLine(std::string line);

public:
    VJAlignmentInfo(std::string filename);

    size_t size() {
        return read_names_.size();
    }

    IgGeneAlignment GetVAlignmentByIndex(size_t index);

    IgGeneAlignment GetJAlignmentByIndex(size_t index);

};