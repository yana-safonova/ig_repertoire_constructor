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

    std::string VJCloneSetDecomposer::GetVKeyByClone(const annotation_utils::AnnotatedClone &clone) const {
        //return clone.VGene().name() + "__" + clone.JGene().name();
        return GetGeneBaseName(clone.VGene().name());
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

    core::Decomposition VJCloneSetDecomposer::CreateDecompositionByVGenes() const {
        std::map<std::string, std::vector<size_t>> v_clusters;
        for(size_t i = 0; i < clone_set_.size(); i++) {
            std::string v_key = GetVKeyByClone(clone_set_[i]);
            std::cout << v_key << "\n";
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

    core::Decomposition VJCloneSetDecomposer::CreateDecompositionToOneClass() const {

        std::map<std::string, std::vector<size_t>> vj_clusters;
        vj_clusters["0"] = std::vector<size_t>();

        for(size_t i = 0; i < clone_set_.size(); i++) {
            vj_clusters["0"].push_back(i);
        }

        core::Decomposition decomposition(clone_set_.size());
        size_t class_index = 0;
        for (auto it = vj_clusters.begin(); it != vj_clusters.end(); it++) {
            auto vj_class = it->second;
            for (auto it2 = vj_class.begin(); it2 != vj_class.end(); it2++)
                decomposition.SetClass(*it2, class_index);
            class_index++;
        }
        return decomposition;
    }


}