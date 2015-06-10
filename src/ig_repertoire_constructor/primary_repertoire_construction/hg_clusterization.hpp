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
        size_t min_shift = min<size_t>(r1.GetFrom(), r2.GetFrom());
        size_t start1 = r1.GetFrom() - min_shift;
        size_t start2 = r2.GetFrom() - min_shift;
        size_t overlap_len = min<size_t>(r1.ReadLength() - start1, r2.ReadLength() - start2);
        for(size_t i = 0; i < overlap_len; i++) {
        	if(r1.FullSequence()[start1 + i] != r2.FullSequence()[start2 + i]) {
                ++dist;
                if (dist > max_distance_)
                    return dist;
            }
        }
        return dist;
    }
};

class HG_DecompositionPrinter {
	HG_DecompositionPtr hg_decomposition_;
	const SplicedReadGroup &read_group_;
	CRS_HammingGraph_Ptr hamming_graph_;
	HG_CollapsedStructs_Ptr collapsed_struct_;

public:
	HG_DecompositionPrinter(HG_DecompositionPtr hg_decomposition,
			const SplicedReadGroup &read_group,
			CRS_HammingGraph_Ptr hamming_graph,
			HG_CollapsedStructs_Ptr collapsed_struct) :
				hg_decomposition_(hg_decomposition),
				read_group_(read_group),
				hamming_graph_(hamming_graph),
				collapsed_struct_(collapsed_struct) { }

	void PrintDecomposition(string output_fname) {
		ofstream out_fhandler(output_fname.c_str());
        //cout << "hamming_graph_->N(): " << hamming_graph_->N() << endl;
        //cout << "Oppa: " << collapsed_struct_->NumberNewVertices() << endl;
        for(size_t i = 0; i < read_group_.size(); i++) {
            //cout << i << endl;
            size_t new_vertex = collapsed_struct_->NewIndexOfOldVertex(i);
            size_t sgraph_id = hg_decomposition_->GetVertexClass(new_vertex);
            out_fhandler << sgraph_id << "\t" << read_group_[i].GetReadName() << "\t" << read_group_[i].GetFrom() << endl;
        }
        out_fhandler.close();
	}
};

template<class HammingDistanceCalculator, class DenseSubgraphConstructor, class DenseSubgraphDecompositor>
class HGClustersConstructor {
    // input data
    size_t max_tau_;
    const ig_config::hg_clusterization_params& params_;

    // auxiliary classes-procedures
    HammingDistanceCalculator calculator_;
    DenseSubgraphConstructor dense_subgraph_constructor_;

    // auxiliary structures
    CRS_HammingGraph_Ptr hamming_graph_ptr_;
    HG_CollapsedStructs_Ptr collapsed_struct_;
    HG_DecompositionPtr dense_sgraphs_decomposition_ptr_;
    HG_DecompositionPtr decomposed_dense_sgraphs_;

    // output struct
    HG_DecompositionPtr final_decomposition_ptr_;

    void InitializeCRSHammingGraph(SplicedReadGroup read_group) {
        TRACE("Initialization of Hamming graph edges");
        vector <HGEdge> hg_edges;
        for (size_t i = 0; i < read_group.size() - 1; i++)
            for (size_t j = i + 1; j < read_group.size(); j++) {
//            	cout << read_group[i].ReadName() << " - " << read_group[j].ReadName() << endl;
                size_t dist = calculator_.HammingDistance(read_group[i], read_group[j]);
//                cout << "Distance: " << dist << endl;
                if (dist <= max_tau_) {
                    //cout << i << " " << j << ", distance: "  << dist << endl;
                    hg_edges.push_back(HGEdge(i, j, dist));
                }
            }
        TRACE("Hamming graph contains " << read_group.size() << " vertices and " << hg_edges.size() << " edges");
        hamming_graph_ptr_ = CRS_HammingGraph_Ptr(new CRS_HammingGraph(read_group.size(), hg_edges));
    }

    void CollapseIdenticalVertices() {
        collapsed_struct_ = HG_CollapsedStructs_Ptr(new HG_CollapsedStructs(hamming_graph_ptr_));
        TRACE("Collapsed graph contains " << collapsed_struct_->NumCollapsedVertices() << " vertices");
    }

    string GetDecompositionFname(size_t group_id) {
        stringstream ss;
        ss << "dense_sgraph_decomposition_" << group_id << "_size_" << hamming_graph_ptr_->N() << ".txt";
        return path::append_path(params_.hgc_io_params.dense_subgraphs_dir, ss.str());
    }

    void OutputDenseSgraphDecomposition(HG_DecompositionPtr dense_sgraph_decomposition, SplicedReadGroup read_group) {
    	if(!params_.hgc_io_params.output_dense_subgraphs)
    		return;
    	string output_fname = GetDecompositionFname(read_group.Id());
    	HG_DecompositionPrinter(dense_sgraph_decomposition, read_group, hamming_graph_ptr_,
    			collapsed_struct_).PrintDecomposition(output_fname);
    	TRACE("Decomposition into dense subgraphs #" << read_group.Id() << " was written to " << output_fname);
    }

    void CreateDenseSubgraphDecomposition(size_t graph_id) {
    	dense_sgraphs_decomposition_ptr_ = dense_subgraph_constructor_.CreateDecomposition(hamming_graph_ptr_,
    			collapsed_struct_, graph_id);
    }

    void DecomposeDenseSubgraphs(SplicedReadGroup read_group) {
        DenseSubgraphDecompositor dense_subgraph_decomposer_(params_.min_recessive_abs_size,
        		params_.min_recessive_rel_size, hamming_graph_ptr_, collapsed_struct_,
        		dense_sgraphs_decomposition_ptr_, read_group);
        decomposed_dense_sgraphs_ = dense_subgraph_decomposer_.CreateDecomposition();
    }

    HG_DecompositionPtr CreateTrivialDecomposition(size_t reads_number) {
        HG_DecompositionPtr decomposition(new HG_Decomposition(reads_number));
        for(size_t i = 0; i < reads_number; i++)
            decomposition->SetClass(i, 0);
        return decomposition;
    }

    bool HammingGraphIsIsolated() {
    	return hamming_graph_ptr_->NZ() == 0 or
    			collapsed_struct_->NumberCollapsedEdges(hamming_graph_ptr_) == 0;
    }

    HG_DecompositionPtr CreateIsolatedDecomposition(size_t reads_number) {
        HG_DecompositionPtr decomposition(new HG_Decomposition(reads_number));
        for(size_t i = 0; i < reads_number; i++)
            decomposition->SetClass(i, i);
        return decomposition;
    }

    void CreateFinalDecomposition() {
    	final_decomposition_ptr_ = HG_DecompositionPtr(new HG_Decomposition(hamming_graph_ptr_->N()));
    	for(size_t i = 0; i < hamming_graph_ptr_->N(); i++) {
    		size_t new_index = collapsed_struct_->NewIndexOfOldVertex(i);
    		size_t new_class = decomposed_dense_sgraphs_->GetVertexClass(new_index);
    		final_decomposition_ptr_->SetClass(i, new_class);
    	}
    }

public:
    HGClustersConstructor(size_t max_tau, const ig_config::hg_clusterization_params& params) :
    	max_tau_(max_tau),
        params_(params),
        calculator_(max_tau),
        dense_subgraph_constructor_(params) { }


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

        if(HammingGraphIsIsolated())
        	return CreateIsolatedDecomposition(read_group.size());

        CreateDenseSubgraphDecomposition(read_group.Id());
        OutputDenseSgraphDecomposition(dense_sgraphs_decomposition_ptr_, read_group);
        TRACE("Decomposition of Hamming graph into dense subgraphs was computed");

        DecomposeDenseSubgraphs(read_group);
        TRACE("Decomposition of dense subgraphs was computed");

        CreateFinalDecomposition();
        TRACE("Final decomposition was computed");

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
