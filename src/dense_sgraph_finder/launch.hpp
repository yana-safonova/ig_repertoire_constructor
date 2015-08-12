#pragma once

#include "graph_utils/sparse_graph.hpp"
#include "graph_utils/graph_collapsed_structure.hpp"
#include "graph_utils/graph_io.hpp"
#include "graph_decomposer/dense_subgraph_constructor.hpp"

namespace dense_subgraph_finder {
    class DenseSubgraphFinder {

    public:
        DenseSubgraphFinder() { }

        int Run(const dsf_config::dense_sgraph_finder_params &dsf_params,
                const dsf_config::io_params &io,
                const dsf_config::metis_io_params &metis_io) {
            INFO("==== Dense subgraph finder starts");
            std::string graph_fname = io.input.graph_filename;
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
                    dsf_params, metis_io, io);
            DecompositionPtr decomposition_ptr = denseSubgraphConstructor.CreateDecomposition(graph_ptr,
                                                                                              collapsed_struct_ptr);
            decomposition_ptr->SaveTo(io.output_base.decomposition_filename);
            INFO("==== Dense subgraph finder ends");
            return 0;
        }
    };
}