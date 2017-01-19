#pragma once

#include "evolutionary_graph_utils/evolutionary_tree.hpp"
#include "evolutionary_annotated_shm.hpp"
#include "tree_based_shm_annotator.hpp"

namespace antevolo {
    class AnnotatedEvolutionaryTree {
        const EvolutionaryTree *tree_ptr_;
        TreeBasedSHMAnnotator shm_annotator_;

        // segment SHMs
        std::map<size_t, std::vector<EvolutionaryAnnotatedSHM> > segment_shms_map_;
        std::vector<EvolutionaryAnnotatedSHM> all_segment_shms_;
        std::set<EvolutionaryAnnotatedSHM, EvolutionaryAnnotatedSHMComparator> unique_segment_shms_;

        // CDR3 SHMs
        std::map<size_t, std::vector<CDR3SHM> > cdr3_shms_map_;
        std::vector<CDR3SHM> all_cdr3_shms_;
        std::set<CDR3SHM, CDR3SHMComparator> unique_cdr3_shms_;

        void CheckConsistencyFatal(size_t clone_id) const;

    public:
        AnnotatedEvolutionaryTree(const EvolutionaryTree &tree) :
                tree_ptr_(&tree),
                shm_annotator_(tree) { }

        void AddSegmentSHMForClone(size_t clone_id, annotation_utils::SHM shm);

        // rel_src_pos and rel_dst_pos - positions wrt to start of CDR3
        void AddCDR3SHMForClone(size_t src_id, size_t dst_id, size_t rel_src_pos, size_t rel_dst_pos);

        const EvolutionaryTree& Tree() const { return *tree_ptr_; }

        size_t NumUniqueSHms() const { return /*all_segment_shms_.size() + all_cdr3_shms_.size();*/
                    unique_segment_shms_.size() + unique_cdr3_shms_.size(); }

        size_t NumSynonymousWrtGermlineSHMs() const;

        size_t NumCDR3SHMs() const { return unique_cdr3_shms_.size(); }

        size_t NumSynonymousSHMs() const { return NumSynonymousSegmentSHMs() + NumSynonymousCDR3SHMs(); }

        size_t NumSynonymousSegmentSHMs() const;

        size_t NumSynonymousCDR3SHMs() const;

        size_t RootDepth() const;

        size_t SHMDepth() const;

        size_t NumAddedSHMs() const; // return numner of SHMs that were added wrt to tree root

        size_t GetNumCDR3SHMsForClone(size_t clone_id) const;

        typedef std::vector<EvolutionaryAnnotatedSHM>::const_iterator SHMConstIterator;

        SHMConstIterator cbegin() const { return all_segment_shms_.cbegin(); }

        SHMConstIterator cend() const { return all_segment_shms_.cend(); }

        size_t GetNumCDR3SHMsFromCloneToRoot(size_t clone_id) const;

        const CloneSetWithFakes& GetCloneSet() const { return tree_ptr_->GetCloneSet(); }
    };
}