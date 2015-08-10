#pragma once

#include "graph_utils/sparse_graph.hpp"
#include "graph_utils/graph_collapsed_structure.hpp"
#include "graph_utils/graph_io.hpp"
#include "graph_decomposer/dense_subgraph_constructor.hpp"

namespace dense_subgraph_finder {
    class DenseSubgraphFinder {

    public:
        DenseSubgraphFinder() { }

        int Run() {
            INFO("==== Dense subgraph finder starts");
            std::string graph_fname = dsf_cfg::get().io.graph_filename;
            GraphReader graph_reader(graph_fname);
            SparseGraphPtr graph_ptr = graph_reader.CreateGraph();
            if (!graph_ptr) {
                INFO("Dense subgraph finder was unable to extract graph from " << graph_fname);
                return 1;
            }
            GraphCollapsedStructurePtr collapsed_struct_ptr = GraphCollapsedStructurePtr(
                    new GraphCollapsedStructure(graph_ptr));
            INFO("Collapsed structure contains " << collapsed_struct_ptr->NumberNewVertices() << " vertices");
            dense_subgraph_finder::MetisDenseSubgraphConstructor denseSubgraphConstructor(
                    dsf_cfg::get().dsf_params, dsf_cfg::get().metis_io, dsf_cfg::get().io.graph_filename);
            DecompositionPtr decomposition_ptr = denseSubgraphConstructor.CreateDecomposition(graph_ptr,
                                                                                              collapsed_struct_ptr);
            decomposition_ptr->SaveTo(dsf_cfg::get().io.decomposition_filename);
            INFO("==== Dense subgraph finder ends");
            return 0;
        }
    };
}