#pragma once

#include "hg_decomposition.hpp"

class GreedyJoiningDecomposition {
    // input parameters
    CRS_HammingGraph_Ptr hamming_graph_ptr_;
    HG_CollapsedStructs_Ptr collapsed_struct_;
    HG_DecompositionPtr basic_decomposition_ptr_;
    double average_fillin_threshold_;

    // auxiliary structs
    map<size_t, set<size_t> > decomposition_graph_;
    vector<bool> class_processed_;
    vector<size_t> class_size_;
    size_t num_processed_;
    vector<size_t> vertex_class_;

    // output parameters
    HG_DecompositionPtr output_decomposition_ptr_;

    void InitializeClassStructs() {
        for(size_t i = 0; i < basic_decomposition_ptr_->Size(); i++)
            class_processed_.push_back(false);
        for(size_t i = 0; i < basic_decomposition_ptr_->Size(); i++)
            class_size_.push_back(basic_decomposition_ptr_->ClassSize(i));
    }

    void InitializeDecompositionGraph() {
        for(size_t i = 0; i < basic_decomposition_ptr_->Size(); i++)
            decomposition_graph_[i] = set<size_t>();
    }

    void InitializeVertexClass() {
        for(size_t i = 0; i < collapsed_struct_->NumberNewVertices(); i++)
            vertex_class_.push_back(size_t(-1));
        for(size_t i = 0; i < basic_decomposition_ptr_->Size(); i++)
            for(auto it = basic_decomposition_ptr_->GetClass(i).begin();
                it != basic_decomposition_ptr_->GetClass(i).end(); it++)
                vertex_class_[*it] = i;
    }

    void Initialize() {
        InitializeClassStructs();
        InitializeDecompositionGraph();
        InitializeVertexClass();
    }

    void UpdateDecompositionGraph(size_t class1, size_t class2) {
        if(class1 == class2)
            return;
        decomposition_graph_[class1].insert(class2);
        decomposition_graph_[class2].insert(class1);
    }

    void CreateDecompositionGraph() {
        for(size_t i = 0; i < hamming_graph_ptr_->N(); i++) {
            size_t v1 = collapsed_struct_->NewIndexOfOldVertex(i);
            for(size_t j = hamming_graph_ptr_->RowIndex()[i]; j < hamming_graph_ptr_->RowIndex()[i + 1]; j++) {
                size_t v2 = collapsed_struct_->NewIndexOfOldVertex(hamming_graph_ptr_->Col()[j]);
                size_t class1 = basic_decomposition_ptr_->GetVertexClass(v1);
                size_t class2 = basic_decomposition_ptr_->GetVertexClass(v2);
                UpdateDecompositionGraph(class1, class2);
            }
        }
    }

    size_t GetMaximalAvailableClass() {
        size_t max_size = 0;
        size_t max_class = size_t(-1);
        for(auto it = decomposition_graph_.begin(); it != decomposition_graph_.end(); it++)
            if(class_size_[it->first] > max_size and !class_processed_[it->first]) {
                max_size = class_size_[it->first];
                max_class = it->first;
            }
        return max_class;
    }

    double ComputeRelativeFillin(size_t class_id, size_t vertex) {
        size_t num_edges_to_class = 0;
        size_t old_vertex = collapsed_struct_->OldVerticesList()[vertex];
        set<size_t> used_main_vertices;

        for(size_t i = hamming_graph_ptr_->RowIndex()[old_vertex];
            i < hamming_graph_ptr_->RowIndex()[old_vertex + 1]; i++) {
            size_t old_neigh = hamming_graph_ptr_->Col()[i];
            size_t new_neigh = collapsed_struct_->NewIndexOfOldVertex(old_neigh);
            size_t neigh_class = vertex_class_[new_neigh];
            if(used_main_vertices.find(new_neigh) == used_main_vertices.end() and
                    neigh_class == class_id) {
                num_edges_to_class++;
                used_main_vertices.insert(new_neigh);
            }
        }
        for(size_t i = hamming_graph_ptr_->RowIndexT()[old_vertex];
            i < hamming_graph_ptr_->RowIndexT()[old_vertex + 1]; i++) {
            size_t old_neigh = hamming_graph_ptr_->ColT()[i];
            size_t new_neigh = collapsed_struct_->NewIndexOfOldVertex(old_neigh);
            size_t neigh_class = vertex_class_[new_neigh];
            if(used_main_vertices.find(new_neigh) == used_main_vertices.end() and
                    neigh_class == class_id) {
                num_edges_to_class++;
                used_main_vertices.insert(new_neigh);
            }
        }
        return double(num_edges_to_class) / double(class_size_[class_id]);
    }

    bool ClassesCanBeGlued(size_t main_class, size_t sec_class) {
        if(class_processed_[sec_class]) {
            TRACE("Class " << sec_class << " is already processed");
            return false;
        }

        auto sec_class_set = basic_decomposition_ptr_->GetClass(sec_class);
        double average_fillin = 0;
        for(auto it = sec_class_set.begin(); it != sec_class_set.end(); it++) {
            double cur_avg_fillin = ComputeRelativeFillin(main_class, *it);
            assert(cur_avg_fillin <= 1);
            average_fillin += cur_avg_fillin;
        }
        average_fillin /= double(sec_class_set.size());
        DEBUG("Primary class: " << main_class << " (" << class_size_[main_class] <<
                "), secondary class: " << sec_class << " (" << class_size_[sec_class] <<
                "), avg edge fillin: " << average_fillin);
        return average_fillin >= average_fillin_threshold_;
    }

    void GlueClasses(size_t main_class, size_t sec_class) {
        // for all neighbours of sec class:
        // (1) add main class
        // (2) remove sec class
        // (3) add neigh in adj list of main class
        for(auto it = decomposition_graph_[sec_class].begin(); it != decomposition_graph_[sec_class].end(); it++) {
            decomposition_graph_[*it].erase(sec_class);
            if(*it != main_class) {
                decomposition_graph_[*it].insert(main_class);
                decomposition_graph_[main_class].insert(*it);
            }
        }
        // erase
        decomposition_graph_.erase(sec_class);
        num_processed_++;

        // todo: add some kind of class id remapping
        class_size_[main_class] += class_size_[sec_class];
        for(auto it = basic_decomposition_ptr_->GetClass(sec_class).begin();
            it != basic_decomposition_ptr_->GetClass(sec_class).end(); it++)
            vertex_class_[*it] = main_class;
    }

    void CreateNewDecomposition() {
        //average_fillin_threshold_ = .8; // todo: compute and move to config
        while(num_processed_ < basic_decomposition_ptr_->Size()) {
            size_t cur_main_class = GetMaximalAvailableClass();
            TRACE("New main class: " << cur_main_class << ", size: " << class_size_[cur_main_class]);
            size_t num_glued = 0;
            auto neigh_set = decomposition_graph_[cur_main_class];
            for(auto it = neigh_set.begin(); it != neigh_set.end(); it++)
                if(ClassesCanBeGlued(cur_main_class, *it))  {
                    GlueClasses(cur_main_class, *it);
                    num_glued++;
                    TRACE("Classes " << cur_main_class << " and " << *it << " were glued");
                }
            if(num_glued == 0) {
                class_processed_[cur_main_class] = true;
                num_processed_++;
            }
            TRACE("Processed " << num_processed_ << " vertices from " << basic_decomposition_ptr_->Size());
            TRACE("-------------");
        }
    }

    void WriteNewDecomposition() {
        // compute max class
        size_t max_class = 0;
        for(size_t i = 0; i < vertex_class_.size(); i++)
            max_class = max<size_t>(max_class, vertex_class_[i]);
        // initialization
        vector<size_t> class_index;
        for(size_t i = 0; i <= max_class; i++)
            class_index.push_back(size_t(-1));
        // mark significant cells
        for(size_t i = 0; i < vertex_class_.size(); i++)
            class_index[vertex_class_[i]] = 0;
        // compute dense index of each class
        size_t index = 0;
        for(size_t i = 0; i <= max_class; i++)
            if(class_index[i] == 0) {
                class_index[i] = index;
                index++;
            }
        // write decomposition in class
        for(size_t i = 0; i < vertex_class_.size(); i++) {
            TRACE("Vertex: " << i << ", class: " << vertex_class_[i] << ", new index of class: " <<
                    class_index[vertex_class_[i]]);
            output_decomposition_ptr_->SetClass(i, class_index[vertex_class_[i]]);
        }
    }

public:
    GreedyJoiningDecomposition(CRS_HammingGraph_Ptr hamming_graph_ptr,
        HG_CollapsedStructs_Ptr collapsed_struct,
        HG_DecompositionPtr basic_decomposition_ptr,
        double average_fillin_threshold) :
            hamming_graph_ptr_(hamming_graph_ptr),
            collapsed_struct_(collapsed_struct),
            basic_decomposition_ptr_(basic_decomposition_ptr),
            average_fillin_threshold_(average_fillin_threshold),
            num_processed_(0),
            output_decomposition_ptr_(new HG_Decomposition(basic_decomposition_ptr_->VertexNumber())) { }

    HG_DecompositionPtr ConstructDecomposition() {
        Initialize();
        DEBUG("Initialization of GreedyJoiningDecomposition");
        CreateDecompositionGraph();
        DEBUG("Decomposition graph was created");
        CreateNewDecomposition();
        DEBUG("New decomposition was computed");
        WriteNewDecomposition();
        DEBUG(output_decomposition_ptr_->Size() << " classes were constructed");
        DEBUG("Maximal class size: " << output_decomposition_ptr_->MaxClassSize());
        return output_decomposition_ptr_;
    }

private:
    DECL_LOGGER("GreedyJoiningDecomposition");
};