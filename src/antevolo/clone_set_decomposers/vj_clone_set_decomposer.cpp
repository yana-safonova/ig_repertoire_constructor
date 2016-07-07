#include "vj_clone_set_decomposer.hpp"
#include <annotation_utils/annotated_clone.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>
#include <vector>

namespace antevolo {
    std::string VJCloneSetDecomposer::GetGeneBaseName(seqan::CharString name) const {
        std::string gene_name = std::string(seqan::toCString(name));
        //return gene_name;
        std::vector<std::string> splits;
        boost::split(splits, gene_name, boost::is_any_of("*"), boost::token_compress_on);
        return splits[0];
    }

    std::string VJCloneSetDecomposer::GetVJKeyByClone(const annotation_utils::AnnotatedClone &clone) const {
        //return clone.VGene().name() + "__" + clone.JGene().name();
        return GetGeneBaseName(clone.VGene().name()) + "_" + GetGeneBaseName(clone.JGene().name());
    }

    core::Decomposition VJCloneSetDecomposer::CreateDecomposition() const {
        std::map<std::string, std::vector<size_t>> vj_clusters;
        for(size_t i = 0; i < clone_set_.size(); i++) {
            std::string vj_key = GetVJKeyByClone(clone_set_[i]);
            if(vj_clusters.find(vj_key) == vj_clusters.end())
                vj_clusters[vj_key] = std::vector<size_t>();
            vj_clusters[vj_key].push_back(i);
        }
        core::Decomposition decomposition(clone_set_.size());
        size_t class_index = 0;
        for(auto it = vj_clusters.begin(); it != vj_clusters.end(); it++) {
            auto vj_class = it->second;
            for(auto it2 = vj_class.begin(); it2 != vj_class.end(); it2++)
                decomposition.SetClass(*it2, class_index);
            class_index++;
        }
        return decomposition;
    }
}