#include "launch.hpp"
#include "graph_utils/graph_splitter.hpp"

namespace {
    class NonParallelDenseSubgraphFinder {
        SparseGraphPtr graph_ptr_;
        GraphCollapsedStructurePtr collapsed_struct_ptr_;
        const dsf_config::dense_sgraph_finder_params &dsf_params_;
        const dsf_config::io_params &io_;
        const dsf_config::metis_io_params &metis_io_;

    public:
        NonParallelDenseSubgraphFinder(SparseGraphPtr graph_ptr,
                                       GraphCollapsedStructurePtr collapsed_struct_ptr,
                                       const dsf_config::dense_sgraph_finder_params &dsf_params,
                                       const dsf_config::io_params &io,
                                       const dsf_config::metis_io_params &metis_io) :
                graph_ptr_(graph_ptr),
                collapsed_struct_ptr_(collapsed_struct_ptr),
                dsf_params_(dsf_params),
                io_(io),
                metis_io_(metis_io) { }

        void Run() {
            dense_subgraph_finder::MetisDenseSubgraphConstructor denseSubgraphConstructor(
                    dsf_params_, metis_io_, io_);
            DecompositionPtr decomposition_ptr = denseSubgraphConstructor.CreateDecomposition(graph_ptr_,
                                                                                              collapsed_struct_ptr_);
            decomposition_ptr->SaveTo(io_.output_base.decomposition_filename);
        }
    };

    class ParallelDenseSubgraphFinder {
        SparseGraphPtr graph_ptr_;
        GraphCollapsedStructurePtr collapsed_struct_ptr_;
        const dsf_config::dense_sgraph_finder_params &dsf_params_;
        const dsf_config::io_params &io_;
        const dsf_config::metis_io_params &metis_io_;

    public:
        ParallelDenseSubgraphFinder(SparseGraphPtr graph_ptr,
                                    GraphCollapsedStructurePtr collapsed_struct_ptr,
                                    const dsf_config::dense_sgraph_finder_params &dsf_params,
                                    const dsf_config::io_params &io,
                                    const dsf_config::metis_io_params &metis_io) :
                graph_ptr_(graph_ptr),
                collapsed_struct_ptr_(collapsed_struct_ptr),
                dsf_params_(dsf_params),
                io_(io),
                metis_io_(metis_io) { }

        void Run() {
            auto connected_components = ConnectedComponentGraphSplitter(graph_ptr_).Split();
            // todo: implement me!
        }
    };
}

int dense_subgraph_finder::DenseSubgraphFinder::Run() {
    INFO("==== Dense subgraph finder starts");
    GraphReader graph_reader(io_.input.graph_filename);
    SparseGraphPtr graph_ptr = graph_reader.CreateGraph();
    if (!graph_ptr) {
        INFO("Dense subgraph finder was unable to extract graph from " << io_.input.graph_filename);
        return 1;
    }
    GraphCollapsedStructurePtr collapsed_struct_ptr = GraphCollapsedStructurePtr(
            new GraphCollapsedStructure(graph_ptr));
    INFO("Collapsed structure contains " << collapsed_struct_ptr->NumberNewVertices() << " vertices");
    if(run_params_.threads_count == 1) {
        INFO("Nonparallel mode was chosen");
        NonParallelDenseSubgraphFinder(graph_ptr, collapsed_struct_ptr, dsf_params_, io_, metis_io_).Run();
    }
    else {
        INFO("Parallel mode was chosen. Number of threads " << run_params_.threads_count);
        ParallelDenseSubgraphFinder(graph_ptr, collapsed_struct_ptr, dsf_params_, io_, metis_io_).Run();
    }
    INFO("==== Dense subgraph finder ends");
    return 0;
}