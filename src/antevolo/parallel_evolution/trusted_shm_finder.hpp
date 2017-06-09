#pragma once

#include "mutation_map.hpp"

namespace antevolo {
    class TrustedSHMFinder {
        const MutationMap &shm_map_;
        const ClonalGraph &clonal_graph_;

        std::map<annotation_utils::SHM, size_t> trusted_shm_mult_;
        std::set<size_t> visited_conflicting_dsts_;
        std::set<size_t> visited_srcs_;

        void AddTrustedSHM(TreeSHM shm);

        void AddAllSHMsFromEdge(std::pair<size_t, size_t> edge);

        void AddSHMsFromConflictingEdges(size_t dst, std::set<size_t> conflicting_src);

        bool SHMIsFound(annotation_utils::SHM shm) const {
            return trusted_shm_mult_.find(shm) != trusted_shm_mult_.end();
        }

        void FatalCheckOfNonzeroTrustedMult();

    public:
        TrustedSHMFinder(const MutationMap &shm_map, const ClonalGraph &clonal_graph) : shm_map_(shm_map),
                                                                                        clonal_graph_(clonal_graph) { }

        void FindTrustedSHMs();


        size_t TrustedSHMNumber() const { return trusted_shm_mult_.size(); }

        size_t TrustedNonTrivialSHMNumber() const;

        size_t TrustedSynonymousNumber() const;

        size_t TrustedHotSpotNumber() const;


        size_t GetTrustedMultiplicity(annotation_utils::SHM shm) const {
            if(SHMIsFound(shm))
                return trusted_shm_mult_.at(shm);
            return 0;
        }

        bool SHMIsTrusted(annotation_utils::SHM shm) const {
            return GetTrustedMultiplicity(shm) > 1;
        }
    };
}