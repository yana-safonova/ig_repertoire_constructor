#include "evolutionary_stats_calculator.hpp"
#include <annotation_utils/shm_comparator.hpp>

namespace antevolo {
    void EvolutionaryStatsCalculator::ComputeTreeStats(EvolutionaryTree tree) {
        size_t total_num_shms = 0;
        for(auto edge = tree.begin(); edge != tree.end(); edge++) {
            auto added_v_shms = annotation_utils::SHMComparator::GetAddedSHMs(clone_set_[edge->src_clone_num].VSHMs(),
                                                                            clone_set_[edge->dst_clone_num].VSHMs());
            auto added_j_shms = annotation_utils::SHMComparator::GetAddedSHMs(clone_set_[edge->src_clone_num].JSHMs(),
                                                                              clone_set_[edge->dst_clone_num].JSHMs());
            total_num_shms += added_v_shms.size();
            total_num_shms += added_j_shms.size();
            // todo: add SHMs in CDR3
        }
        size_t root_id = tree.GetRoot();
        total_num_shms += clone_set_[root_id].VSHMs().size();
        total_num_shms += clone_set_[root_id].JSHMs().size();
        INFO("Total numner of SHMs in tree: " << total_num_shms);
    }
}