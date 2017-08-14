#include "cluster_fillin_calculator.hpp"

namespace antevolo {

    size_t ClusterFillinCalculator::CalculateE(boost::unordered_set<size_t> cluster) {
        auto edge_constructor = GetEdgeConstructor();
        size_t non_orphans_number = 0;
        for (size_t i : cluster) {
            bool i_is_orphan = true;
            const auto& clone_i = clone_set_[i];
            for (size_t j : cluster) {
                if (i == j) {
                    continue;
                }
                auto edge = edge_constructor->ConstructEdge(clone_set_[j], clone_i, j, i);
                if (i_is_orphan && edge->IsDirected()) {
                    i_is_orphan = false;
                    ++non_orphans_number;
                }
            }
        }
        return non_orphans_number;
    }
}