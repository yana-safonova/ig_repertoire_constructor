#pragma once

#include "../include_me.hpp"
#include "../spliced_read.hpp"

#include "hamming_graph_clusterization/crs_matrix.hpp"
#include "hamming_graph_clusterization/permutation.hpp"
#include "hamming_graph_clusterization/hg_decomposition.hpp"

#include "dense_subgraph_constructor.hpp"
#include "dense_subgraph_decomposer.hpp"

#include "primary_repertoire_constructor.hpp"
#include "../../include/path_helper.hpp"

namespace ig_repertoire_constructor {

class SimpleHammingDistanceCalculator {
    size_t max_distance_;
public:
    SimpleHammingDistanceCalculator(size_t max_distance) :
            max_distance_(max_distance) {
    }

    size_t HammingDistance(const SplicedRead &r1, const SplicedRead &r2) {
        size_t dist = 0;
        for (size_t i = 0; i < r1.size(); ++i)
            if (r1[i] != r2[i]) {
                ++dist;
                if (dist > max_distance_)
                    return dist;
            }
        return dist;
    }
};

template<class HammingDistanceCalculator, class DenseSubgraphConstructor, class DenseSubgraphDecompositor>
class HGClustersConstructor {
    // auxiliary classes-procedures
    HammingDistanceCalculator calculator_;
    DenseSubgraphConstructor dense_subgraph_constructor_;
    DenseSubgraphDecompositor dense_subgraph_decomposer_;

    // input data
    size_t max_tau_;

    // auxiliary structures
    CRS_HammingGraph_Ptr hamming_graph_ptr_;
    HG_CollapsedStructs_Ptr collapsed_struct_;
    HG_DecompositionPtr dense_sgraphs_decomposition_ptr_;

    // output struct
    HG_DecompositionPtr final_decomposition_ptr_;

    void InitializeCRSHammingGraph(SplicedReadGroup read_group) {
        TRACE("Initialization of Hamming graph edges");
        vector <HGEdge> hg_edges;
        for (size_t i = 0; i < read_group.size() - 1; i++)
            for (size_t j = i + 1; j < read_group.size(); j++) {
                size_t dist = calculator_.HammingDistance(read_group[i], read_group[j]);
                if (dist <= max_tau_)
                    hg_edges.push_back(HGEdge(i, j, dist));
            }
        TRACE("Hamming graph contains " << read_group.size() << " vertices and " << hg_edges.size() << " edges");
        hamming_graph_ptr_ = CRS_HammingGraph_Ptr(new CRS_HammingGraph(hg_edges));
    }

    void CollapseIdenticalVertices() {
        collapsed_struct_ = HG_CollapsedStructs_Ptr(new HG_CollapsedStructs(hamming_graph_ptr_));
        TRACE("Collapsed graph contains " << collapsed_struct_->NumCollapsedVertices() << " vertices");
    }

    void CreateDenseSubgraphDecomposition(size_t graph_id) {
    	dense_sgraphs_decomposition_ptr_ = dense_subgraph_constructor_.CreateDecomposition(hamming_graph_ptr_,
    			collapsed_struct_, graph_id);
    }

    void DecomposeDenseSubgraphs(SplicedReadGroup read_group) {
    	final_decomposition_ptr_ = dense_subgraph_decomposer_.CreateDecomposition(hamming_graph_ptr_,
    			collapsed_struct_, dense_sgraphs_decomposition_ptr_, read_group);
    }

    /* void ComputeHGFinalDecomposition(SplicedReadGroup read_group) {
        TRACE("AlignmentDecompositionConstructor starts");
        AlignmentDecompositionConstructor align_constructor(hamming_graph_ptr_, collapsed_struct_,
                dense_sgraphs_decomposition_ptr_, read_group);
        final_decomposition_ptr_ = align_constructor.ConstructDecomposition();

        TRACE("DecompositionStatsCalculator starts");
        DecompositionStatsCalculator calculator(final_decomposition_ptr_, hamming_graph_ptr_, collapsed_struct_);
        calculator.WriteStatsInFile("final_decomposition_stats.txt");
    }

    string GetDecompositionFname(size_t group_id, string prefix) {
        stringstream ss;
        ss << prefix << "_decomposition_" << group_id;
        if(prefix == "secondary")
            ss << "_" << class_joining_edge_threshold_;
        ss << ".txt";
        return path::append_path(ig_cfg::get().io.hgraph_dir, ss.str());
    }

    void SaveDecomposition(const vector <SplicedRead> &spliced_reads, HG_DecompositionPtr decomposition,
            string filename) {
        ofstream out(filename.c_str());
        for(size_t i = 0; i < spliced_reads.size(); i++) {
            size_t old_vertex = i;
            size_t new_vertex = collapsed_struct_->NewIndexOfOldVertex(i);
            size_t sgraph_id = decomposition->GetVertexClass(new_vertex);
            out << old_vertex << "\t" << sgraph_id << "\t" << spliced_reads[i].GetReadName() << "\t" <<
                    spliced_reads[i].GetFrom() << endl;
        }
        out.close();
    }
    */

    HG_DecompositionPtr CreateTrivialDecomposition(size_t reads_number) {
        HG_DecompositionPtr decomposition(new HG_Decomposition(reads_number));
        for(size_t i = 0; i < reads_number; i++)
            decomposition->SetClass(i, 0);
        return decomposition;
    }

public:
    HGClustersConstructor(size_t max_tau,
    		double edge_perc_threshold,
    		double class_joining_edge_threshold,
    		size_t min_recessive_abs_size,
    		double min_recessive_rel_size) :
        calculator_(max_tau),
        dense_subgraph_constructor_(edge_perc_threshold, class_joining_edge_threshold),
        dense_subgraph_decomposer_(min_recessive_abs_size, min_recessive_rel_size),
        max_tau_(max_tau) { }


    HG_DecompositionPtr ConstructClusters(SplicedReadGroup read_group) {
        if(read_group.size() < 4) {
            return CreateTrivialDecomposition(read_group.size());
        }

        TRACE("------------");
        TRACE("Clusterization starts. Index " << read_group.Id());

        InitializeCRSHammingGraph(read_group);
        TRACE("CRS Hamming graph was constructed");

        CollapseIdenticalVertices();
        TRACE("Identical vertices were collapsed");

        CreateDenseSubgraphDecomposition(read_group.Id());
        TRACE("Decomposition of Hamming graph into dense subgraphs was computed");

        DecomposeDenseSubgraphs(read_group);
        TRACE("Decomposition of dense subgraphs was computed");

        return final_decomposition_ptr_;
    }

    HG_CollapsedStructs_Ptr CollapsedVerticesStruct() {
        return collapsed_struct_;
    }

private:
    DECL_LOGGER("HGClustersConstructor");
};

typedef HGClustersConstructor<SimpleHammingDistanceCalculator, MetisDenseSubgraphConstructor, IterativeDenseSubgraphDecomposer> StandardHGClustersConstructor;

}
