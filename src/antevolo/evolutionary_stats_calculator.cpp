#include "evolutionary_stats_calculator.hpp"
#include <annotation_utils/shm_comparator.hpp>

namespace antevolo {
    AnnotatedEvolutionaryTree EvolutionaryStatsCalculator::ComputeTreeStats(const EvolutionaryTree &tree) {
        AnnotatedEvolutionaryTree annotated_tree(clone_set_, tree);
        for(auto edge = tree.cbegin(); edge != tree.cend(); edge++) {
            auto added_v_shms = annotation_utils::SHMComparator::GetAddedSHMs(clone_set_[edge->src_clone_num].VSHMs(),
                                                                            clone_set_[edge->dst_clone_num].VSHMs());
            auto added_j_shms = annotation_utils::SHMComparator::GetAddedSHMs(clone_set_[edge->src_clone_num].JSHMs(),
                                                                              clone_set_[edge->dst_clone_num].JSHMs());
            for(auto it = added_v_shms.begin(); it != added_v_shms.end(); it++)
                annotated_tree.AddSHMForClone(edge->dst_clone_num, *it);
            for(auto it = added_j_shms.begin(); it != added_j_shms.end(); it++)
                annotated_tree.AddSHMForClone(edge->dst_clone_num, *it);
            // todo: add SHMs in CDR3
        }
//        size_t root_id = tree.GetRoot();
//        total_num_shms += clone_set_[root_id].VSHMs().size();
//        total_num_shms += clone_set_[root_id].JSHMs().size();
        INFO("Total number of SHMs in a tree with " << tree.NumVertives() << " vertices: " <<
                     annotated_tree.NumUniqueSHms());
        return annotated_tree;
    }
}