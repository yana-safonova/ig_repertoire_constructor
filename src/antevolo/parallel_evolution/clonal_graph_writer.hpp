#pragma once

#include "clonal_graph.hpp"

namespace antevolo {
    class ClonalGraphWriter {
        const ClonalGraph &cgraph_;

        //std::string GetEdgeLabel(const ParallelRhomb &rhomb, RhombSide rhomb_side, size_t index) const;

    public:
        ClonalGraphWriter(const ClonalGraph &cgraph) : cgraph_(cgraph) {
        }

        void operator()(std::string output_fname);
    };
}