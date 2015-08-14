#include <omp.h>
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
                    dsf_params_,
                    metis_io_,
                    io_.output_nonparallel.graph_copy_filename,
                    io_.output_base.decomposition_filename);
            DecompositionPtr decomposition_ptr = denseSubgraphConstructor.CreateDecomposition(graph_ptr_,
                                                                                              collapsed_struct_ptr_);
            INFO("Dense subgraph decomposition was written to " << io_.output_base.decomposition_filename);
        }
    };

    class ParallelDenseSubgraphFinder {
        SparseGraphPtr graph_ptr_;
        GraphCollapsedStructurePtr collapsed_struct_ptr_;
        const dsf_config::dense_sgraph_finder_params &dsf_params_;
        const dsf_config::io_params &io_;
        const dsf_config::metis_io_params &metis_io_;

        string GetSubgraphFilename(size_t subgraph_index) {
            stringstream ss;
            ss << "subgraph_" << subgraph_index << ".graph";
            return path::append_path(io_.output_mthreading.connected_components_dir, ss.str());
        }

        string GetDecompositionFilename(size_t subgraph_index) {
            stringstream ss;
            ss << "decomposition_" << subgraph_index << ".txt";
            return path::append_path(io_.output_mthreading.decompositions_dir, ss.str());
        }

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
            vector<SparseGraphPtr> connected_components = ConnectedComponentGraphSplitter(graph_ptr_).Split();
#pragma omp parallel for
            for(size_t i = 0; i < connected_components.size(); i++) {
                SparseGraphPtr current_subgraph = connected_components[i];
                GraphCollapsedStructurePtr current_collapsed_struct = GraphCollapsedStructurePtr(
                        new GraphCollapsedStructure(current_subgraph));
                string graph_filename = GetSubgraphFilename(i);
                string decomposition_filename = GetDecompositionFilename(i);
                dense_subgraph_finder::MetisDenseSubgraphConstructor denseSubgraphConstructor(
                        dsf_params_,
                        metis_io_,
                        graph_filename,
                        decomposition_filename);
                DecompositionPtr decomposition_ptr = denseSubgraphConstructor.CreateDecomposition(graph_ptr_,
                                                                                                  collapsed_struct_ptr_);
                INFO("Dense subgraph decomposition was written to " << decomposition_filename);
            }
            INFO("Parallel construction of dense subgraphs for connected components finished");
            INFO("Connected components in GRAPH format were written to " <<
                         io_.output_mthreading.connected_components_dir);
            INFO("Dense subgraph decomposition for connected components were written to " <<
                         io_.output_mthreading.decompositions_dir);
        }
        // todo: create final decomposition
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
        omp_set_num_threads(run_params_.threads_count);
        NonParallelDenseSubgraphFinder(graph_ptr, collapsed_struct_ptr, dsf_params_, io_, metis_io_).Run();
    }
    else {
        INFO("Parallel mode was chosen. Number of threads: " << run_params_.threads_count);
        ParallelDenseSubgraphFinder(graph_ptr, collapsed_struct_ptr, dsf_params_, io_, metis_io_).Run();
    }
    INFO("==== Dense subgraph finder ends");
    return 0;
}