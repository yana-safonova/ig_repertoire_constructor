#pragma once

#include "../graph_utils/sparse_graph.hpp"
#include <boost/unordered_set.hpp>

namespace antevolo {

    class CDR3HammingGraphComponentInfo {
        typedef std::map<std::string, std::vector<size_t>> UniqueCDR3IndexMap;
        typedef std::map<std::string, size_t> CDR3ToIndexMap;

        GraphComponentMap& graph_component_map_;
        const UniqueCDR3IndexMap& cdr3_to_indices_vector_map_;
        const CDR3ToIndexMap& cdr3_to_old_index_map_;
        const std::vector<std::string>& unique_cdr3s_;
        SparseGraphPtr hg_component_;
        size_t component_id_;

    public:
        CDR3HammingGraphComponentInfo(GraphComponentMap& graph_component_map,
                             const UniqueCDR3IndexMap& cdr3_to_indices_vector_map,
                             const CDR3ToIndexMap& cdr3_to_old_index_map,
                             const std::vector<std::string>& unique_cdr3s,
                             SparseGraphPtr hg_component,
                             size_t component_id) :
            graph_component_map_(graph_component_map),
            cdr3_to_indices_vector_map_(cdr3_to_indices_vector_map),
            cdr3_to_old_index_map_(cdr3_to_old_index_map),
            unique_cdr3s_(unique_cdr3s),
            hg_component_(hg_component),
            component_id_(component_id) {}

        size_t GetOldIndexByNewIndex(size_t i);
        size_t GetNewIndexByOldIndex(size_t old_index);
        const std::vector<size_t>& GetClonesByCDR3(std::string cdr3);
        const std::string& GetCDR3ByOldIndex(size_t old_index);
        size_t GetOldIndexByCDR3(const std::string& cdr3);
        const std::vector<size_t>& GetClonesByOldIndex(size_t old_index);
        SparseGraph::EdgesIterator GetSimilarCDR3sBeginByOldIndex(size_t old_index);
        SparseGraph::EdgesIterator GetSimilarCDR3sEndByOldIndex(size_t old_index);
        boost::unordered_set<size_t> GetAllClones();
    };

}

