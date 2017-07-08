#pragma once

#include "mutation_map.hpp"
#include "trusted_shm_finder.hpp"

namespace antevolo {
    class MutationMapWriter {
        const MutationMap &map_;
        //const TrustedSHMFinder &trusted_finder_;

        std::string EdgesToString(std::vector<std::pair<size_t, size_t> > edges) const;

        std::string TripletPairToString(TreeSHM shm) const;

    public:
        MutationMapWriter(const MutationMap &map) : map_(map) { }

        void operator()(std::string output_fname) const;
    };
}