#pragma once

#include "evolutionary_graph_utils/evolutionary_tree.hpp"
#include "../shm_counting/tree_shm_map.hpp"
#include "../shm_counting/tree_shm_calculator.hpp"
#include "../clone_set_with_fakes.hpp"

namespace antevolo {
    class AnnotatedEvolutionaryTree {
        const CloneSetWithFakes &clone_set_;
        const EvolutionaryTree &tree_;

        // statistics
        UniqueSHMCalculator shm_calculator_;
        TreeSHMMap shm_map_;
        std::map<size_t, size_t> clone_added_shm_map_; // clone id -> number of added SHMs WRT parent clone
        std::map<size_t, size_t> clone_added_from_root_map_; // clone id -> number of added SHMs WRT root

        void CheckConsistencyFatal(size_t clone_id) const;

        void InitializeCloneSHMMap();

    public:
        AnnotatedEvolutionaryTree(const EvolutionaryTree &tree) : clone_set_(tree.GetCloneSet()),
                                                                  tree_(tree),
                                                                  shm_calculator_(clone_set_, tree),
                                                                  shm_map_( ) { //clone_set_, tree_){
            for(auto it = tree_.cbegin(); it != tree_.cend(); it++) {
                shm_calculator_.AddSHMsFromEdge(*it);
            }
            shm_map_  = shm_calculator_.GetSHMMap();
            InitializeCloneSHMMap();
        }

        const EvolutionaryTree& Tree() const { return tree_; }

        const TreeSHMMap &SHMMap() const { return shm_map_; }


        size_t NumUniqueSHMs() const { return shm_map_.size(); }

        size_t NumSynonymousSHMs() const { return shm_map_.NumSynonymousSHMs(); }

        size_t NumSHMsInCDRs() const { return shm_map_.NumCDRSHMs(); }

        size_t NumSHMsInRegion(annotation_utils::StructuralRegion region) const {
            return shm_map_.NumSHMsInRegion(region);
        }

        size_t NumSynonymousSHMsInRegion(annotation_utils::StructuralRegion region) const {
            return shm_map_.NumSynSHMsInRegion(region);
        }

        size_t MaxSHMMultiplicity() const { return shm_map_.MaxMultiplicity(); }

        size_t RootDepth() const;

        size_t TreeDepth() const;

        size_t GetRegionLength(annotation_utils::StructuralRegion region) const;
    };
}