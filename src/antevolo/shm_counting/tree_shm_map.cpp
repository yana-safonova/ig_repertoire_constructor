#include "tree_shm_map.hpp"

namespace antevolo {
    void TreeSHMMap::AddSHM(TreeSHM shm) {
        if(shm_mult_map_.find(shm) == shm_mult_map_.end())
            shm_mult_map_[shm] = 0;
        shm_mult_map_[shm]++;
    }
}