#include "tree_shm_map.hpp"

namespace antevolo {
    void TreeSHMMap::AddSHM(TreeSHM shm, size_t src_id, size_t dst_id) {
        if(shm_mult_map_.find(shm) == shm_mult_map_.end())
            shm_mult_map_[shm] = 0;
        shm_mult_map_[shm]++;
        if(shm_clone_ids_.find(shm) == shm_clone_ids_.end())
            shm_clone_ids_[shm] = std::vector<std::pair<size_t, size_t> >();
        shm_clone_ids_[shm].push_back(std::make_pair(src_id, dst_id));
    }

    size_t TreeSHMMap::NumSynonymousSHMs() const {
        size_t num_synonymous = 0;
        for(auto it = shm_mult_map_.cbegin(); it != shm_mult_map_.cend(); it++) {
            if(it->first.Synonymous())
                num_synonymous++;
        }
        return num_synonymous;
    }

    size_t TreeSHMMap::NumCDRSHMs() const {
        return NumSHMsInRegion(annotation_utils::StructuralRegion::CDR1) +
                NumSHMsInRegion(annotation_utils::StructuralRegion::CDR2) +
                NumSHMsInRegion(annotation_utils::StructuralRegion::CDR3);
    }

    size_t TreeSHMMap::MaxMultiplicity() const {
        size_t max_mult = 0;
        for(auto it = shm_mult_map_.cbegin(); it != shm_mult_map_.cend(); it++)
            max_mult = std::max(max_mult, it->second);
        return max_mult;
    }

    size_t TreeSHMMap::NumSHMsInRegion(annotation_utils::StructuralRegion region) const {
        size_t num_region_shms = 0;
        for(auto it = shm_mult_map_.cbegin(); it != shm_mult_map_.cend(); it++) {
            if(it->first.region == region)
                num_region_shms++;
        }
        return num_region_shms;
    }

    size_t TreeSHMMap::NumSynSHMsInRegion(annotation_utils::StructuralRegion region) const {
        size_t num_syn_shms = 0;
        for(auto it = shm_mult_map_.cbegin(); it != shm_mult_map_.cend(); it++) {
            if(it->first.region == region and it->first.Synonymous())
                num_syn_shms++;
        }
        return num_syn_shms;
    }
}