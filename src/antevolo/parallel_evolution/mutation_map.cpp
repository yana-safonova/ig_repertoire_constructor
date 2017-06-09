#include "mutation_map.hpp"

namespace antevolo {
    void MutationMap::FillMutationFromEdge(size_t src, size_t dst) {
        auto added_shms = annotation_utils::SHMComparator::GetAddedSHMs(clone_set_[src].VSHMs(),
                                                                        clone_set_[dst].VSHMs());
        std::vector<TreeSHM> subs_shms;
        TreeSHMCalculator tree_shm_calc(clone_set_);
        for(auto it = added_shms.begin(); it != added_shms.end(); it++) {
            if(it->shm_type != annotation_utils::SHMType::SubstitutionSHM)
                continue;
            subs_shms.push_back(tree_shm_calc.ComputeTreeSHM(*it, src, dst));
        }
        edge_shms_map_[std::make_pair(src, dst)] = subs_shms;
    }

    void MutationMap::FillMutationsFromGraph() {
        for(auto it = clonal_graph_.cbegin(); it != clonal_graph_.cend(); it++) {
            auto src = it->first;
            auto dst_set = it->second;
            for(auto v = dst_set.begin(); v != dst_set.end(); v++) {
                FillMutationFromEdge(src, *v);
            }
        }
    }

    void MutationMap::FillSHMEdgeMap() {
        for(auto it = edge_shms_map_.begin(); it != edge_shms_map_.end(); it++) {
            auto edge = it->first;
            auto shms = it->second;
            for(auto s = shms.begin(); s != shms.end(); s++) {
                if(shm_edge_map_.find(*s) == shm_edge_map_.end())
                    shm_edge_map_[*s] = std::vector<std::pair<size_t, size_t> >();
                shm_edge_map_[*s].push_back(edge);
            }
        }
    }

    std::set<TreeSHM> MutationMap::GetUniqueOutgoingSHMs(size_t src) const {
        auto outgoing_v = clonal_graph_.OutgoingVertices(src);
        std::set<TreeSHM> unique_shms;
        for(auto v = outgoing_v.begin(); v != outgoing_v.end(); v++) {
            auto edge_shms = edge_shms_map_.at(std::make_pair(src, *v));
            for(auto s = edge_shms.begin(); s != edge_shms.end(); s++) {
                unique_shms.insert(*s);
            }
        }
        return unique_shms;
    }

    void MutationMap::ComputeUniqueMutationMult() {
        for(auto it = clonal_graph_.cbegin(); it != clonal_graph_.cend(); it++) {
            auto src = it->first;
            auto unique_outgoing_shms = GetUniqueOutgoingSHMs(src);
            for(auto s = unique_outgoing_shms.begin(); s != unique_outgoing_shms.end(); s++) {
                if(unique_shm_mult_.find(*s) == unique_shm_mult_.end())
                    unique_shm_mult_[*s] = 0;
                unique_shm_mult_[*s]++;
            }
        }
    }

    size_t MutationMap::NumNontrivialSHMs() const {
        size_t num_nontrivial = 0;
        for(auto it = unique_shm_mult_.begin(); it != unique_shm_mult_.end(); it++) {
            if(it->second > 1)
                num_nontrivial++;
        }
        return num_nontrivial;
    }

    size_t MutationMap::NumSynonymousNonTrivialSHMs() const {
        size_t num_synonymous = 0;
        for(auto it = unique_shm_mult_.begin(); it != unique_shm_mult_.end(); it++) {
            if(it->second == 1)
                continue;
            if(it->first.Synonymous())
                num_synonymous++;
        }
        return num_synonymous;
    }

    size_t MutationMap::MaxMultiplicity() const {
        size_t max_multiplicity = 0;
        for (auto it = unique_shm_mult_.begin(); it != unique_shm_mult_.end(); it++) {
            max_multiplicity = std::max(max_multiplicity, it->second);
        }
        return max_multiplicity;
    }

    bool MutationMap::SHMIsHotspot(TreeSHM shm, std::pair<size_t, size_t> edge) const {
        auto read_seq = clone_set_[edge.second].Read().seq;
        // R = A/G Y = C/T W = A/T
        bool rgyw = read_seq[shm.src_pos - 1] == 'A' or read_seq[shm.src_pos - 1] == 'G';
        rgyw = rgyw and (read_seq[shm.src_pos + 1] == 'C' or read_seq[shm.src_pos + 1] == 'T');
        rgyw = rgyw and (read_seq[shm.src_pos + 2] == 'A' or read_seq[shm.src_pos + 2] == 'T');
        bool wrcy = read_seq[shm.src_pos - 2] == 'A' or read_seq[shm.src_pos - 2] == 'T';
        wrcy = wrcy and (read_seq[shm.src_pos - 1] == 'A' or read_seq[shm.src_pos - 1] == 'G');
        wrcy = wrcy and (read_seq[shm.src_pos + 1] == 'C' or read_seq[shm.src_pos + 1] == 'T');
        return rgyw or wrcy;
    }

    size_t MutationMap::NumNonTrivialHotSpotSHMs() const {
        size_t num_hotspots = 0;
        for(auto it = shm_edge_map_.begin(); it != shm_edge_map_.end(); it++) {
            if(!SHMIsNonTrivial(it->first))
                continue;
            auto first_edge = it->second[0];
            if(SHMIsHotspot(it->first, first_edge))
                num_hotspots++;
        }
        return num_hotspots;
    }

    bool MutationMap::SHMIsNonTrivial(TreeSHM shm) const {
        VERIFY_MSG(unique_shm_mult_.find(shm) != unique_shm_mult_.end(), "SHM " << shm << " was not found");
        return unique_shm_mult_.at(shm) > 1;
    }

    bool MutationMap::SHMIsHotSpot(TreeSHM shm) const {
        VERIFY(shm_edge_map_.find(shm) != shm_edge_map_.end());
        auto edge = shm_edge_map_.at(shm)[0];
        return SHMIsHotspot(shm, edge);
    }

    size_t MutationMap::GetSHMMultiplicity(TreeSHM shm) const {
        VERIFY(unique_shm_mult_.find(shm) != unique_shm_mult_.end());
        return unique_shm_mult_.at(shm);
    }

    std::vector<TreeSHM > MutationMap::GetSHMsByEdge(size_t src, size_t dst) const {
        auto edge_pair = std::make_pair(src, dst);
        VERIFY(edge_shms_map_.find(edge_pair) != edge_shms_map_.end());
        return edge_shms_map_.at(edge_pair);
    }

    std::vector<std::pair<size_t, size_t> > MutationMap::GetEdgesBySHM(TreeSHM shm) const {
        VERIFY(shm_edge_map_.find(shm) != shm_edge_map_.end());
        return shm_edge_map_.at(shm);
    }

    seqan::CharString MutationMap::VGeneName() const {
        for(auto it = edge_shms_map_.cbegin(); it != edge_shms_map_.cend(); it++)
            return clone_set_[it->first.first].VAlignment().subject().name();
    }
}