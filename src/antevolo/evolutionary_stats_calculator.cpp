#include "evolutionary_stats_calculator.hpp"
#include <annotation_utils/shm_comparator.hpp>

namespace antevolo {
    AnnotatedEvolutionaryTree EvolutionaryStatsCalculator::AddSHMsFromRoot(AnnotatedEvolutionaryTree annotated_tree) const {
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

    AnnotatedEvolutionaryTree EvolutionaryStatsCalculator::AddSHMsFromEdges(AnnotatedEvolutionaryTree annotated_tree) const {
        const EvolutionaryTree& tree = annotated_tree.Tree();
        for (auto edge = tree.cbegin(); edge != tree.cend(); edge++) {
            //std::cout << edge->dst_clone_num << " -> " << edge->src_clone_num << std::endl;
            auto added_v_shms = annotation_utils::SHMComparator::GetAddedSHMs(clone_set_[edge->src_clone_num].VSHMs(),
                                                                              clone_set_[edge->dst_clone_num].VSHMs());
            auto added_j_shms = annotation_utils::SHMComparator::GetAddedSHMs(clone_set_[edge->src_clone_num].JSHMs(),
                                                                              clone_set_[edge->dst_clone_num].JSHMs());
            //std::cout << added_v_shms.size() << " - " << added_j_shms.size() << std::endl;
            for (auto it = added_v_shms.begin(); it != added_v_shms.end(); it++)
                annotated_tree.AddSHMForClone(edge->dst_clone_num, *it);
            for (auto it = added_j_shms.begin(); it != added_j_shms.end(); it++)
                annotated_tree.AddSHMForClone(edge->dst_clone_num, *it);
            //std::cout << clone_set_[edge->src_clone_num].CDR3() << " - " <<
            //        clone_set_[edge->dst_clone_num].CDR3() << std::endl;
            // todo: add SHMs in CDR3
        }
        return annotated_tree;
    }

    AnnotatedEvolutionaryTree EvolutionaryStatsCalculator::ComputeStatsForTree(const EvolutionaryTree &tree) const {
        AnnotatedEvolutionaryTree annotated_tree(clone_set_, tree);
        if(tree.IsForest()) {
            // todo: after implementation of tree splitting turn this into assert
            INFO("Input tree is a forest containing " << tree.GetRootNumber() <<
                         " trees. Computation of stats will be skipped since preliminary decomposition is required");
            return annotated_tree;
        }
        annotated_tree = AddSHMsFromRoot(annotated_tree);
        annotated_tree = AddSHMsFromEdges(annotated_tree);
        INFO("# vertices: " << tree.NumEdges() + 1 << ", root depth: " << annotated_tree.RootDepth() <<
                     ", # unique SHMs: " << annotated_tree.NumUniqueSHms() <<
                     ", # added SHMs: " << annotated_tree.NumAddedSHMs() <<
                     ", # synonymous SHMs: " << annotated_tree.NumSynonymousSHMs() <<
                     ", # synonymous SHMs wrt to germline: " << annotated_tree.NumSynonymousWrtGermlineSHMs());
        return annotated_tree;
    }

    AnnotatedTreeStorage EvolutionaryStatsCalculator::ComputeStatsForStorage(
            const EvolutionaryTreeStorage &storage) const {
        AnnotatedTreeStorage annotated_storage(clone_set_);
        for(auto it = storage.cbegin(); it != storage.cend(); it++) {
            annotated_storage.AddAnnotatedTree(ComputeStatsForTree(*it));
        }
        return annotated_storage;
    }
}