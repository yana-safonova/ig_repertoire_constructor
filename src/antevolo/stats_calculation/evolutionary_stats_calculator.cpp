#include "evolutionary_stats_calculator.hpp"
#include <annotation_utils/shm_comparator.hpp>
#include "evolutionary_graph_utils/evolutionary_tree_splitter.hpp"

namespace antevolo {
    AnnotatedEvolutionaryTree EvolutionaryStatsCalculator::AddSHMsFromRoot(AnnotatedEvolutionaryTree annotated_tree) const {
        const auto& clone_set = annotated_tree.GetCloneSet();
        auto root_id = annotated_tree.Tree().GetRoot();
        auto v_shms = clone_set[root_id].VSHMs();
        for(auto it = v_shms.cbegin(); it != v_shms.cend(); it++) {
            annotated_tree.AddSegmentSHMForClone(root_id, *it);
        }
        auto j_shms = clone_set[root_id].JSHMs();
        for(auto it = j_shms.cbegin(); it != j_shms.cend(); it++) {
            annotated_tree.AddSegmentSHMForClone(root_id, *it);
        }
        return annotated_tree;
    }

    AnnotatedEvolutionaryTree EvolutionaryStatsCalculator::AddSHMsFromEdges(AnnotatedEvolutionaryTree annotated_tree) const {
        const auto& clone_set = annotated_tree.GetCloneSet();
        const EvolutionaryTree& tree = annotated_tree.Tree();
        for (auto it = tree.cbegin(); it != tree.cend(); it++) {
            auto edge = *it;
            auto added_v_shms = annotation_utils::SHMComparator::GetAddedSHMs(clone_set[edge->SrcNum()].VSHMs(),
                                                                              clone_set[edge->DstNum()].VSHMs());
            auto added_j_shms = annotation_utils::SHMComparator::GetAddedSHMs(clone_set[edge->SrcNum()].JSHMs(),
                                                                              clone_set[edge->DstNum()].JSHMs());
            for (auto it = added_v_shms.begin(); it != added_v_shms.end(); it++)
                annotated_tree.AddSegmentSHMForClone(edge->DstNum(), *it);
            for (auto it = added_j_shms.begin(); it != added_j_shms.end(); it++)
                annotated_tree.AddSegmentSHMForClone(edge->DstNum(), *it);
            if(clone_set[edge->SrcNum()].CDR3() != clone_set[edge->DstNum()].CDR3()) {
                auto src_cdr3 = clone_set[edge->SrcNum()].CDR3();
                auto dst_cdr3 = clone_set[edge->DstNum()].CDR3();
                // todo: works only for mismatches in CDR3, refactor this in future
                VERIFY_MSG(seqan::length(src_cdr3) == seqan::length(dst_cdr3), "CDR3 " << src_cdr3 <<
                        " and " << dst_cdr3 << " have different lengths");
                for(size_t i = 0; i < seqan::length(src_cdr3); i++) {
                    if(src_cdr3[i] != dst_cdr3[i]) {
                        //std::cout << clone_set[edge->src_clone_num].CDR3Range().start_pos << std::endl;
                        //std::cout << clone_set[edge->dst_clone_num].CDR3Range().start_pos << std::endl;
                        annotated_tree.AddCDR3SHMForClone(edge->SrcNum(), edge->DstNum(), i, i);
                    }
                }
            }
        }
        return annotated_tree;
    }

    AnnotatedEvolutionaryTree EvolutionaryStatsCalculator::ComputeStatsForTree(const EvolutionaryTree &tree) const {
        AnnotatedEvolutionaryTree annotated_tree(tree);
        VERIFY_MSG(!tree.IsForest(), "Input tree is a forest containing " << tree.GetRootNumber() << " trees");
        annotated_tree = AddSHMsFromRoot(annotated_tree);
        annotated_tree = AddSHMsFromEdges(annotated_tree);
        //INFO("# vertices: " << tree.NumEdges() + 1 << ", root depth: " << annotated_tree.RootDepth() <<
        //     ", # unique SHMs: " << annotated_tree.NumUniqueSHms() <<
        //     ", # added SHMs: " << annotated_tree.NumAddedSHMs() <<
        //     ", # synonymous SHMs: " << annotated_tree.NumSynonymousSHMs() <<
        //     ", # synonymous SHMs wrt to germline: " << annotated_tree.NumSynonymousWrtGermlineSHMs());
        return annotated_tree;
    }

    AnnotatedTreeStorage EvolutionaryStatsCalculator::ComputeStatsForStorage(
            const EvolutionaryTreeStorage &storage) const {
        AnnotatedTreeStorage annotated_storage;
        for(auto it = storage.cbegin(); it != storage.cend(); it++) {
            annotated_storage.AddAnnotatedTree(ComputeStatsForTree(*it));
        }
        return annotated_storage;
    }
}