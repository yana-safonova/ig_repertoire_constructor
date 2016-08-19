#include "evolutionary_stats_calculator.hpp"
#include <annotation_utils/shm_comparator.hpp>

namespace antevolo {
    AnnotatedEvolutionaryTree EvolutionaryStatsCalculator::AddSHMsFromRoot(AnnotatedEvolutionaryTree annotated_tree) {
        auto root_id = annotated_tree.Tree().GetRoot();
        auto v_shms = clone_set_[root_id].VSHMs();
        for(auto it = v_shms.cbegin(); it != v_shms.cend(); it++) {
            annotated_tree.AddSHMForClone(root_id, *it);
        }
        auto j_shms = clone_set_[root_id].JSHMs();
        for(auto it = j_shms.cbegin(); it != j_shms.cend(); it++) {
            annotated_tree.AddSHMForClone(root_id, *it);
        }
        return annotated_tree;
    }

    AnnotatedEvolutionaryTree EvolutionaryStatsCalculator::AddSHMsFromEdges(AnnotatedEvolutionaryTree annotated_tree) {
        const EvolutionaryTree& tree = annotated_tree.Tree();
        for (auto edge = tree.cbegin(); edge != tree.cend(); edge++) {
            auto added_v_shms = annotation_utils::SHMComparator::GetAddedSHMs(clone_set_[edge->src_clone_num].VSHMs(),
                                                                              clone_set_[edge->dst_clone_num].VSHMs());
            auto added_j_shms = annotation_utils::SHMComparator::GetAddedSHMs(clone_set_[edge->src_clone_num].JSHMs(),
                                                                              clone_set_[edge->dst_clone_num].JSHMs());
            for (auto it = added_v_shms.begin(); it != added_v_shms.end(); it++)
                annotated_tree.AddSHMForClone(edge->dst_clone_num, *it);
            for (auto it = added_j_shms.begin(); it != added_j_shms.end(); it++)
                annotated_tree.AddSHMForClone(edge->dst_clone_num, *it);
            // todo: add SHMs in CDR3
        }
        return annotated_tree;
    }

    AnnotatedEvolutionaryTree EvolutionaryStatsCalculator::ComputeTreeStats(const EvolutionaryTree &tree) {
        AnnotatedEvolutionaryTree annotated_tree(clone_set_, tree);
        if(tree.IsForest()) {
            // todo: after implementation of tree splitting turn this into assert
            INFO("Input tree is a forest containing " << tree.GetRootNumber() <<
                         " trees. Computation of stats will be skipped since preliminary decomposition is required");
            return annotated_tree;
        }
        annotated_tree = AddSHMsFromRoot(annotated_tree);
        annotated_tree = AddSHMsFromEdges(annotated_tree);
        INFO("Total number of SHMs in a tree with " << tree.NumVertives() << " vertices: " <<
                     annotated_tree.NumUniqueSHms() << ". # synonymous SHMs: " <<
                     annotated_tree.NumSynonymousSHMs());
        return annotated_tree;
    }
}