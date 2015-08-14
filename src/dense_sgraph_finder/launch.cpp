#include <omp.h>
#include "launch.hpp"
#include "graph_utils/graph_splitter.hpp"

namespace {
    class NonParallelDenseSubgraphFinder {
        SparseGraphPtr graph_ptr_;
        const dsf_config::dense_sgraph_finder_params &dsf_params_;
        const dsf_config::io_params &io_;
        const dsf_config::metis_io_params &metis_io_;

    public:
        NonParallelDenseSubgraphFinder(SparseGraphPtr graph_ptr,
                                       const dsf_config::dense_sgraph_finder_params &dsf_params,
                                       const dsf_config::io_params &io,
                                       const dsf_config::metis_io_params &metis_io) :
                graph_ptr_(graph_ptr),
                dsf_params_(dsf_params),
                io_(io),
                metis_io_(metis_io) { }

        void Run() {
            dense_subgraph_finder::MetisDenseSubgraphConstructor denseSubgraphConstructor(
                    dsf_params_,
                    metis_io_,
                    io_.output_nonparallel.graph_copy_filename,
                    io_.output_base.decomposition_filename);
            GraphCollapsedStructurePtr collapsed_struct_ptr = GraphCollapsedStructurePtr(
                    new GraphCollapsedStructure(graph_ptr_));
            cout << *collapsed_struct_ptr << endl;
            INFO("Collapsed structure contains " << collapsed_struct_ptr->NumberNewVertices() << " vertices");

            DecompositionPtr decomposition_ptr = denseSubgraphConstructor.CreateDecomposition(graph_ptr_,
                                                                                              collapsed_struct_ptr);
            INFO("Dense subgraph decomposition was written to " << io_.output_base.decomposition_filename);
        }
    };

    class ParallelDenseSubgraphFinder {
        SparseGraphPtr graph_ptr_;
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

        // this function is incorrect
        // todo: fix me!!!
        DecompositionPtr CreateFinalDecomposition(size_t num_connected_components) {
            GraphComponentMap &component_map = graph_ptr_->GetGraphComponentMap();
            cout << component_map << endl;
            GraphCollapsedStructurePtr collapsed_structure(new GraphCollapsedStructure(graph_ptr_));
            map<size_t, size_t> vertex_new_set;
            size_t cur_set_id = 0;
            for(size_t i = 0; i < num_connected_components; i++) {
                cout << "1: " << i << endl;
                string cur_decomposition_fname = GetDecompositionFilename(i);
                cout << "2: " << cur_decomposition_fname << endl;
                DecompositionPtr subgraph_decomposition(new Decomposition(cur_decomposition_fname));
                cout << "3: " << subgraph_decomposition->Size() << endl;
                for(size_t j = 0; j < subgraph_decomposition->Size(); j++) {
                    const set<size_t> &cur_subclass = subgraph_decomposition->GetClass(j);
                    for(auto it = cur_subclass.begin(); it != cur_subclass.end(); it++) {
                        cout << "4: " << j << " " << *it << endl;
                        size_t new_vertex = *it;
                        size_t old_vertex = component_map.GetOldVertexByNewVertex(i, new_vertex);
                        cout << "5: " << old_vertex << endl;
                        size_t new_set_id = cur_set_id + j;
                        vertex_new_set[old_vertex] = new_set_id;
                    }
                }
                cur_set_id += subgraph_decomposition->Size();
            }
            DecompositionPtr final_decomposition_ptr(new Decomposition(graph_ptr_->N()));
            for(size_t i = 0; i < graph_ptr_->N(); i++) {
                size_t new_vertex_id = collapsed_structure->NewIndexOfOldVertex(i);
                assert(vertex_new_set.find(new_vertex_id) != vertex_new_set.end());
                size_t class_id = vertex_new_set[new_vertex_id];
                final_decomposition_ptr->SetClass(i, class_id);
            }
            cout << "Final decomposition: " << endl;
            cout << *final_decomposition_ptr;
            return final_decomposition_ptr;
        }

    public:
        ParallelDenseSubgraphFinder(SparseGraphPtr graph_ptr,
                                    const dsf_config::dense_sgraph_finder_params &dsf_params,
                                    const dsf_config::io_params &io,
                                    const dsf_config::metis_io_params &metis_io) :
                graph_ptr_(graph_ptr),
                dsf_params_(dsf_params),
                io_(io),
                metis_io_(metis_io) { }

        void Run() {
            vector<SparseGraphPtr> connected_components = ConnectedComponentGraphSplitter(graph_ptr_).Split();
//#pragma omp parallel for
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
                DecompositionPtr decomposition_ptr = denseSubgraphConstructor.CreateDecomposition(current_subgraph,
                                                                                                  current_collapsed_struct);
                INFO("Dense subgraph decomposition was written to " << decomposition_filename);
            }
            //assert(false);
            INFO("Parallel construction of dense subgraphs for connected components finished");
            INFO("Connected components in GRAPH format were written to " <<
                         io_.output_mthreading.connected_components_dir);
            INFO("Dense subgraph decomposition for connected components were written to " <<
                         io_.output_mthreading.decompositions_dir);
            DecompositionPtr final_decomposition = CreateFinalDecomposition(connected_components.size());
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
    if(run_params_.threads_count == 1) {
        INFO("Nonparallel mode was chosen");
        NonParallelDenseSubgraphFinder(graph_ptr, dsf_params_, io_, metis_io_).Run();
    }
    else {
        INFO("Parallel mode was chosen. Number of threads: " << run_params_.threads_count);
        omp_set_num_threads(run_params_.threads_count);
        ParallelDenseSubgraphFinder(graph_ptr, dsf_params_, io_, metis_io_).Run();
    }
    INFO("==== Dense subgraph finder ends");
    return 0;
}