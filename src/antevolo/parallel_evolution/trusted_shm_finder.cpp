#include "trusted_shm_finder.hpp"

namespace antevolo {

    void TrustedSHMFinder::AddTrustedSHM(annotation_utils::SHM shm) {
        if(trusted_shm_mult_.find(shm) == trusted_shm_mult_.end())
            trusted_shm_mult_[shm] = 0;
        trusted_shm_mult_[shm]++;
    }

    void TrustedSHMFinder::AddAllSHMsFromEdge(std::pair<size_t, size_t> edge) {
        if(visited_srcs_.find(edge.first) != visited_srcs_.end())
            return;
        auto unique_outgoing_shms = shm_map_.GetUniqueOutgoingSHMs(edge.first);
        for(auto s = unique_outgoing_shms.begin(); s != unique_outgoing_shms.end(); s++) {
            AddTrustedSHM(*s);
        }
        visited_srcs_.insert(edge.first);
    }

    void TrustedSHMFinder::AddSHMsFromConflictingEdges(size_t dst, std::set<size_t> conflicting_src) {
        if(visited_conflicting_dsts_.find(dst) != visited_conflicting_dsts_.end())
            return;
        std::map<annotation_utils::SHM, size_t> shm_mult;
        for(auto v = conflicting_src.begin(); v != conflicting_src.end(); v++) {
            auto shms = shm_map_.GetSHMsByEdge(*v, dst);
            for(auto s = shms.begin(); s != shms.end(); s++) {
                if(shm_mult.find(*s) == shm_mult.end())
                    shm_mult[*s] = 0;
                shm_mult[*s]++;
            }
        }
        for(auto it = shm_mult.begin(); it != shm_mult.end(); it++)
            if(it->second == conflicting_src.size())
                AddTrustedSHM(it->first);
        visited_conflicting_dsts_.insert(dst);
    }

    void TrustedSHMFinder::FatalCheckOfNonzeroTrustedMult() {
        for(auto it = shm_map_.shm_mult_cbegin(); it != shm_map_.shm_mult_cend(); it++) {
            if(GetTrustedMultiplicity(it->first) != 0)
                continue;
            auto edges = shm_map_.GetEdgesBySHM(it->first);
            std::cout << "== ERROR: " << it->first << std::endl;
            for(auto e = edges.begin(); e != edges.end(); e++) {
                std::cout << e->first << " - " << e->second << "; ";
            }
            std::cout << std::endl;
        }
    }

    void TrustedSHMFinder::FindTrustedSHMs() {
        for(auto it = shm_map_.shm_mult_cbegin(); it != shm_map_.shm_mult_cend(); it++) {
            auto edges = shm_map_.GetEdgesBySHM(it->first);
            for(auto e = edges.begin(); e != edges.end(); e++) {
               if(!clonal_graph_.EdgeIsEndOfBulge(e->first, e->second)) {
                   AddAllSHMsFromEdge(*e);
               }
               else {
                   auto conflicting_edges = clonal_graph_.GetConflictingEdges(e->second);
                   VERIFY(conflicting_edges.size() > 1);
                   AddSHMsFromConflictingEdges(e->second, conflicting_edges);
               }
            }
        }
        FatalCheckOfNonzeroTrustedMult();
    }

    size_t TrustedSHMFinder::TrustedNonTrivialSHMNumber() const {
        size_t non_trivial_num = 0;
        for(auto it = trusted_shm_mult_.cbegin(); it != trusted_shm_mult_.cend(); it++) {
            if(!SHMIsTrusted(it->first))
                continue;
            if(shm_map_.SHMIsNonTrivial(it->first))
                non_trivial_num++;
        }
        return non_trivial_num;
    }

    size_t TrustedSHMFinder::TrustedSynonymousNumber() const {
        size_t synonymous_num = 0;
        for(auto it = trusted_shm_mult_.cbegin(); it != trusted_shm_mult_.cend(); it++) {
            if(!SHMIsTrusted(it->first))
                continue;
            if(shm_map_.SHMIsNonTrivial(it->first) and shm_map_.SHMIsSynonymous(it->first))
                synonymous_num++;
        }
        return synonymous_num;
    }

    size_t TrustedSHMFinder::TrustedHotSpotNumber() const {
        size_t hotspot_num = 0;
        for(auto it = trusted_shm_mult_.cbegin(); it != trusted_shm_mult_.cend(); it++) {
            if(!SHMIsTrusted(it->first))
                continue;
            if(shm_map_.SHMIsNonTrivial(it->first) and shm_map_.SHMIsHotSpot(it->first))
                hotspot_num++;
        }
        return hotspot_num;
    }
}