#pragma once

#include "hg_decomposition.hpp"

class DecompositionStatsCalculator {
    HG_DecompositionPtr decomposition_;
    CRS_HammingGraph_Ptr hamming_graph_;
    HG_CollapsedStructs_Ptr collapsed_struct_;

    vector<size_t> num_edges_in_class;

    size_t GetClassOfVertices(size_t v1, size_t v2) {
        if(decomposition_->GetVertexClass(v1) == decomposition_->GetVertexClass(v2))
            return decomposition_->GetVertexClass(v1);
        return size_t(-1);
    }

    void Initialize() {
        for(size_t i = 0; i < decomposition_->Size(); i++)
            num_edges_in_class.push_back(0);
        for(size_t i = 0; i < hamming_graph_->N(); i++)
            for(size_t j = hamming_graph_->RowIndex()[i]; j < hamming_graph_->RowIndex()[i + 1]; j++) {
                size_t v1 = i;
                size_t v2 = hamming_graph_->Col()[j];
                if(!collapsed_struct_->VertexIsMain(v1) or !collapsed_struct_->VertexIsMain(v2))
                    continue;
                // if both vertices are main
                size_t v1_new = collapsed_struct_->NewVerticesList()[v1];
                size_t v2_new = collapsed_struct_->NewVerticesList()[v2];
                size_t class_id = GetClassOfVertices(v1_new, v2_new);
                if(class_id == size_t(-1))
                    continue;
                num_edges_in_class[class_id]++;
            }
    }

    double ComputeClassFillin(size_t class_id) {
        if(decomposition_->ClassSize(class_id) == 1)
            return 0;
        return double(num_edges_in_class[class_id]) / double(decomposition_->ClassSize(class_id) *
                (decomposition_->ClassSize(class_id) - 1)) * 2;
    }

    size_t ComputeClassSize(size_t class_id) {
        return decomposition_->RealSizeOfClass(class_id, collapsed_struct_);

        size_t class_size = 0;
        for(auto it = decomposition_->GetClass(class_id).begin();
            it != decomposition_->GetClass(class_id).end(); it++) {
            size_t old_index = collapsed_struct_->OldVerticesList()[*it];
            class_size += collapsed_struct_->GetMultiplicityOf(old_index);
        }
        return class_size;
    }

public:
    DecompositionStatsCalculator(HG_DecompositionPtr decomposition,
            CRS_HammingGraph_Ptr hamming_graph,
            HG_CollapsedStructs_Ptr collapsed_struct) :
            decomposition_(decomposition),
            hamming_graph_(hamming_graph),
            collapsed_struct_(collapsed_struct) { }

    void WriteStatsInFile(string filename) {
        TRACE("Initialization starts");
        Initialize();
        TRACE("Initialization ends");
        ofstream out(filename.c_str(), std::ios_base::app);
        for(size_t i = 0; i < decomposition_->Size(); i++) {
            assert(decomposition_->ClassSize(i) != 0);
            out << i << "\t" << decomposition_->ClassSize(i) << "\t" <<
                    ComputeClassSize(i) << "\t" << ComputeClassFillin(i) << endl;
        }
        TRACE("Statistics were appended to " << filename);
        out.close();
    }

private:
    DECL_LOGGER("DecompositionStatsCalculator");
};