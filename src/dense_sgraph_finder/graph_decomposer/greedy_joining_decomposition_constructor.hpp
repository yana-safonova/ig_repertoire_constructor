#pragma once

#include "../graph_utils/decomposition.hpp"
#include "../graph_utils/sparse_graph.hpp"

namespace dense_subgraph_finder {

    class GreedyJoiningDecomposition {
        // input parameters
        SparseGraphPtr hamming_graph_ptr_;
        DecompositionPtr basic_decomposition_ptr_;
        double average_fillin_threshold_;
        size_t min_supernode_size_;

        // auxiliary structs
        map <size_t, set<size_t>> decomposition_graph_;
        vector <bool> class_processed_;
        vector <bool> class_has_supernode_;
        vector <size_t> class_size_;
        size_t num_processed_;
        vector <size_t> vertex_class_;

        // output parameters
        DecompositionPtr output_decomposition_ptr_;

        void InitializeClassStructs();

        void InitializeDecompositionGraph();

        void InitializeVertexClass();

        void Initialize();

        void UpdateDecompositionGraph(size_t class1, size_t class2);

        void CreateDecompositionGraph();

        size_t GetMaximalAvailableClass();

        double ComputeRelativeFillin(size_t class_id, size_t vertex);

        bool ClassesCanBeGlued(size_t main_class, size_t sec_class);

        void GlueClasses(size_t main_class, size_t sec_class);

        void CreateNewDecomposition();

        void WriteNewDecomposition();

    public:
        GreedyJoiningDecomposition(SparseGraphPtr hamming_graph_ptr,
                                   DecompositionPtr basic_decomposition_ptr,
                                   double average_fillin_threshold,
                                   size_t min_supernode_size) :
                hamming_graph_ptr_(hamming_graph_ptr),
                basic_decomposition_ptr_(basic_decomposition_ptr),
                average_fillin_threshold_(average_fillin_threshold),
                min_supernode_size_(min_supernode_size),
                num_processed_(0),
                output_decomposition_ptr_(new Decomposition(basic_decomposition_ptr_->VertexNumber())) { }

        DecompositionPtr ConstructDecomposition();

    private:
        DECL_LOGGER("GreedyJoiningDecomposition");
    };

}