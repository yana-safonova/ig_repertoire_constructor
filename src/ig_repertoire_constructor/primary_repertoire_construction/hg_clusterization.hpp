#pragma once

#include "../include_me.hpp"
#include "../spliced_read.hpp"

#include "hamming_graph_clusterization/crs_matrix.hpp"
#include "hamming_graph_clusterization/permutation.hpp"
#include "hamming_graph_clusterization/hg_decomposition.hpp"

#include "hamming_graph_clusterization/simple_decomposition_constructor.hpp"
#include "hamming_graph_clusterization/alignment_decomposition_constructor.hpp"
#include "hamming_graph_clusterization/greedy_joining_decomposition_constructor.hpp"
#include "hamming_graph_clusterization/decomposition_stats_calculator.hpp"

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

template<class HammingDistanceCalculator>
class HGClustersConstructor {
    // auxiliary classes-procedures
    HammingDistanceCalculator calculator_;

    // input data
    size_t max_tau_;
    double edge_perc_threshold_;
    double class_joining_edge_threshold_;

    // temporary data
    CRS_HammingGraph_Ptr hamming_graph_ptr_;

    // auxiliary structures
    HG_CollapsedStructs_Ptr collapsed_struct_;
    PermutationPtr permutation_ptr_;

    HG_DecompositionPtr primary_decomposition_ptr_;
    HG_DecompositionPtr secondary_decomposition_ptr_;
    HG_DecompositionPtr final_decomposition_ptr_;

    void InitializeCRSHammingGraph(const vector <SplicedRead> &spliced_reads) {
        TRACE("Initialization of Hamming graph edges");
        vector <HGEdge> hg_edges;
        for (size_t i = 0; i < spliced_reads.size() - 1; i++)
            for (size_t j = i + 1; j < spliced_reads.size(); j++) {
                size_t dist = calculator_.HammingDistance(spliced_reads[i], spliced_reads[j]);
                if (dist <= max_tau_)
                    hg_edges.push_back(HGEdge(i, j, dist));
            }
//        for (auto it = hg_edges.begin(); it != hg_edges.end(); it++)
//            TRACE(it->i << " " << it->j << " " << it->dist);
        TRACE("Hamming graph contains " << spliced_reads.size() << " vertices and " << hg_edges.size() << " edges");
        hamming_graph_ptr_ = CRS_HammingGraph_Ptr(new CRS_HammingGraph(hg_edges));
    }

    void CollapseIdenticalVertices() {
        collapsed_struct_ = HG_CollapsedStructs_Ptr(new HG_CollapsedStructs(hamming_graph_ptr_));
        TRACE("Collapsed graph contains " << collapsed_struct_->NumCollapsedVertices() << " vertices");
    }

    string GetPermutationFname(size_t group_id) {
        stringstream ss;
        ss << "hamming_graphs_8_collapsed_tau3/hgraph_" << group_id << "_size_" <<
                hamming_graph_ptr_->N() << "_nsize_" <<
                collapsed_struct_->NumCollapsedVertices() << "_tau_3.graph.iperm";
        return ss.str();
    }

    bool PermutationExists(size_t group_id) {
        string perm_fname = GetPermutationFname(group_id);
        // temporary
        return !ifstream(perm_fname).fail();
    }

    void InitializePermutation(size_t group_id) {
        string perm_fname = GetPermutationFname(group_id);
        TRACE("Permutation will be extracted from " << perm_fname);
        permutation_ptr_ = PermutationPtr(new Permutation(collapsed_struct_->NumCollapsedVertices()));
        permutation_ptr_->ReadFromFile(perm_fname);
    }

    void ComputeHGPrimaryDecomposition() {
        SimpleDecompositionConstructor simple_constructor(hamming_graph_ptr_, permutation_ptr_, collapsed_struct_,
                edge_perc_threshold_);
        primary_decomposition_ptr_ = simple_constructor.CreateDecomposition();

        // temporary
        DecompositionStatsCalculator calculator(primary_decomposition_ptr_, hamming_graph_ptr_, collapsed_struct_);
        calculator.WriteStatsInFile("primary_decomposition_stats.txt");
    }

    void ImprovePrimaryDecomposition() {
        TRACE("GreedyJoiningDecomposition threshold: " << class_joining_edge_threshold_);
        GreedyJoiningDecomposition decomposition_improver(hamming_graph_ptr_, collapsed_struct_,
                primary_decomposition_ptr_, class_joining_edge_threshold_);
        secondary_decomposition_ptr_ = decomposition_improver.ConstructDecomposition();

        // temporary
        TRACE("DecompositionStatsCalculator starts");
        DecompositionStatsCalculator calculator(secondary_decomposition_ptr_, hamming_graph_ptr_, collapsed_struct_);
        stringstream ss;
        ss << "secondary_decomposition_stats_" << class_joining_edge_threshold_ << ".txt";
        calculator.WriteStatsInFile(ss.str());
    }

    void ComputeHGFinalDecomposition(const vector <SplicedRead> &spliced_reads, size_t group_id) {
        TRACE("AlignmentDecompositionConstructor starts");
        AlignmentDecompositionConstructor align_constructor(hamming_graph_ptr_, collapsed_struct_,
                secondary_decomposition_ptr_, spliced_reads, group_id);
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

    HG_DecompositionPtr CreateTrivialDecomposition(size_t reads_number) {
        HG_DecompositionPtr decomposition(new HG_Decomposition(reads_number));
        for(size_t i = 0; i < reads_number; i++)
            decomposition->SetClass(i, 0);
        return decomposition;
    }

public:
    HGClustersConstructor(size_t max_tau, double edge_perc_threshold, double class_joining_edge_threshold) :
        calculator_(max_tau),
        max_tau_(max_tau),
        edge_perc_threshold_(edge_perc_threshold),
        class_joining_edge_threshold_(class_joining_edge_threshold) { }


    HG_DecompositionPtr ConstructClusters(const vector <SplicedRead> &spliced_reads, size_t group_id) {
        if(spliced_reads.size() < 4) {
            return CreateTrivialDecomposition(spliced_reads.size());
        }

        TRACE("------------");
        TRACE("Clusterization starts. Index " << group_id);

        InitializeCRSHammingGraph(spliced_reads);
        TRACE("CRS Hamming graph was constructed");

        CollapseIdenticalVertices();
        TRACE("Identical vertices were collapsed");

        if(!PermutationExists(group_id)) {
            TRACE("Graph is too small or too big. Permutation was not found");
            return CreateTrivialDecomposition(spliced_reads.size());
        }
        InitializePermutation(group_id);
        TRACE("Permutation was initialized");
        //TRACE(*permutation_ptr_);

        ComputeHGPrimaryDecomposition();
        TRACE("Primary HG decomposition was computed");
        //TRACE(*primary_decomposition_ptr_);
        SaveDecomposition(spliced_reads, primary_decomposition_ptr_, GetDecompositionFname(group_id, "primary"));
        //primary_decomposition_ptr_->SaveTo(GetDecompositionFname(group_id, "primary"));
        TRACE("Primary decomposition was saved to " << GetDecompositionFname(group_id, "primary"));

        TRACE("Improvement of primary decomposition");
        ImprovePrimaryDecomposition();
        TRACE("Secondary decomposition was computed");
        //TRACE(*secondary_decomposition_ptr_);
        SaveDecomposition(spliced_reads, secondary_decomposition_ptr_, GetDecompositionFname(group_id, "secondary"));
        TRACE("Primary decomposition was saved to " << GetDecompositionFname(group_id, "secondary"));

        ComputeHGFinalDecomposition(spliced_reads, group_id);
        TRACE("Final HG decomposition was computed");
        //TRACE(*final_decomposition_ptr_);
        SaveDecomposition(spliced_reads, final_decomposition_ptr_, GetDecompositionFname(group_id, "final"));
        TRACE("Primary decomposition was saved to " << GetDecompositionFname(group_id, "final"));

        return final_decomposition_ptr_;
    }

    HG_CollapsedStructs_Ptr CollapsedVerticesStruct() {
        return collapsed_struct_;
    }

private:
    DECL_LOGGER("HGClustersConstructor");
};

}
