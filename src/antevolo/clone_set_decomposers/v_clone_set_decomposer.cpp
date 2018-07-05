#include "v_clone_set_decomposer.hpp"

#include "vj_clone_set_decomposer.hpp"
#include <annotation_utils/annotated_clone.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>
#include <vector>

namespace antevolo {


    std::string VCloneSetDecomposer::GetGeneKeyByClone(const annotation_utils::AnnotatedClone& clone) const {
        return GetGeneBaseName(clone.VGene().name());
    }

    core::Decomposition VCloneSetDecomposer::CreateDecomposition() const {
        std::map<std::string, std::vector<size_t>> v_clusters;
        for(size_t i = 0; i < clone_set_.size(); i++) {
            std::string v_key = GetGeneKeyByClone(clone_set_[i]);
            if(v_clusters.find(v_key) == v_clusters.end())
                v_clusters[v_key] = std::vector<size_t>();
            v_clusters[v_key].push_back(i);
        }
        core::Decomposition decomposition(clone_set_.size());
        size_t class_index = 0;
        for(auto it = v_clusters.begin(); it != v_clusters.end(); it++) {
            auto v_class = it->second;
            for(auto it2 = v_class.begin(); it2 != v_class.end(); it2++)
                decomposition.SetClass(*it2, class_index);
            class_index++;
        }
        return decomposition;
    }
}