#pragma once

#include "permutation.hpp"
#include "crs_matrix.hpp"

class HG_Decomposition {
    // input params
    size_t num_vertices_;

    // decomposition fields
    vector<size_t> vertex_class_;
    vector<set<size_t> > decomposition_classes_;

    void InitializeVertexClasses() {
        for(size_t i = 0; i < num_vertices_; i++)
            vertex_class_.push_back(0);
    }

    void AddNewClass() {
        set<size_t> new_class;
        decomposition_classes_.push_back(new_class);
    }

public:
    HG_Decomposition(size_t num_vertices) :
            num_vertices_(num_vertices) {
        InitializeVertexClasses();
    }

    void SetClass(size_t vertex, size_t class_id) {
        assert(vertex < num_vertices_);
        if(decomposition_classes_.size() == 0 or class_id > decomposition_classes_.size() - 1)
            for(size_t i = decomposition_classes_.size(); i <= class_id; i++)
                AddNewClass();

        decomposition_classes_[class_id].insert(vertex);
        vertex_class_[vertex] = class_id;
    }

    const set<size_t>& GetClass(size_t index) const {
        assert(index < Size());
        return decomposition_classes_[index];
    }

    const set<size_t>& LastClass() const {
        return GetClass(Size() - 1);
    }

    size_t ClassSize(size_t index) const {
        return GetClass(index).size();
    }

    size_t LastClassSize() const {
        return LastClass().size();
    }

    size_t VertexNumber() const {
        return vertex_class_.size();
    }

    size_t GetVertexClass(size_t vertex) const {
        assert(vertex < VertexNumber());
        return vertex_class_[vertex];
    }

    size_t Size() const { return decomposition_classes_.size(); }

    void SaveTo(string output_fname) {
        TRACE("Decomposition will written to " << output_fname);
        ofstream out(output_fname.c_str());
        for(auto it = vertex_class_.begin(); it != vertex_class_.end(); it++)
            out << *it << endl;
        out.close();
    }

    bool LastClassContains(size_t vertex) const {
        return LastClass().find(vertex) != LastClass().end();
    }

    size_t MaxClassSize() {
        size_t max_class = 0;
        for(auto it = decomposition_classes_.begin(); it != decomposition_classes_.end(); it++)
            max_class = max<size_t>(max_class, it->size());
        return max_class;
    }

    bool IsTrivial() { return Size() == 1; }

    size_t RealSizeOfClass(size_t class_id, HG_CollapsedStructs_Ptr collapsed_struct_ptr) {
        assert(class_id < Size());
        size_t class_size = 0;
        for(auto it = decomposition_classes_[class_id].begin();
            it != decomposition_classes_[class_id].end(); it++) {
            size_t old_index = collapsed_struct_ptr->OldVerticesList()[*it];
            class_size += collapsed_struct_ptr->GetMultiplicityOf(old_index);
        }
        return class_size;
    }

private:
    DECL_LOGGER("HG_Decomposition");
};

ostream& operator<<(ostream &out, const HG_Decomposition &hg_decomposition) {
    out << "Decomposition contains " << hg_decomposition.Size() << " classes" << endl;
    for(size_t i = 0; i < hg_decomposition.Size(); i++) {
        out << "Class " << i << ", size: " << hg_decomposition.ClassSize(i) << ". ";
        for(auto it = hg_decomposition.GetClass(i).begin(); it != hg_decomposition.GetClass(i).end(); it++)
            out << *it << " ";
        out << endl;
    }
    return out;
}

typedef shared_ptr<HG_Decomposition> HG_DecompositionPtr;

/*
Computaton of statistics
    size_t GetEdgeClass(size_t v1, size_t v2) {
        for(size_t i = 0; i < decomposition_sets_.size(); i++) {
            auto cur_set = decomposition_sets_[i];
            if(cur_set.find(v1) != cur_set.end() and cur_set.find(v2) != cur_set.end())
                return i;
            if(cur_set.find(v1) != cur_set.end() or cur_set.find(v2) != cur_set.end())
                break;
        }
        return size_t(-1);
    }

    vector<size_t> ClassesFilling() {
        vector<size_t> num_edges_inside_class;
        for(size_t i = 0; i < decomposition_sets_.size(); i++)
            num_edges_inside_class.push_back(0);

        //for(size_t i = 0; i < collapsed_struct_->MainVerticesMap().size(); i++)
        //    cout << i << " " << collapsed_struct_->MainVerticesMap()[i] << " ";
        //cout << endl;

        for(size_t i = 0; i < hamming_graph_ptr_->N(); i++)
            for(size_t j = hamming_graph_ptr_->RowIndex()[i]; j < hamming_graph_ptr_->RowIndex()[i + 1]; j++) {
                size_t row_old = i;
                size_t col_old = hamming_graph_ptr_->Col()[j];
                if(!collapsed_struct_->VertexIsMain(row_old) or !collapsed_struct_->VertexIsMain(col_old))
                    continue;
                size_t row_new = collapsed_struct_->NewVerticesList()[collapsed_struct_->MainVerticesMap()[row_old]];
                size_t col_new = collapsed_struct_->NewVerticesList()[collapsed_struct_->MainVerticesMap()[col_old]];
                assert(collapsed_struct_->NewIndexValid(row_new) and collapsed_struct_->NewIndexValid(col_new));
                size_t class_id = GetEdgeClass(row_new, col_new);
                if(class_id != size_t(-1))
                    num_edges_inside_class[class_id]++;
           }
        return num_edges_inside_class;
    }

    double RelativeFilling(const vector<size_t>& classes_fillin) {
        size_t num_edges_inside_classes = 0;
        for(auto it = classes_fillin.begin(); it != classes_fillin.end(); it++)
            num_edges_inside_classes += *it;
        return double(num_edges_inside_classes) / double(hamming_graph_ptr_->NZ());
    }

    vector<size_t> ComputeRealClassSizes() {
        vector <size_t> real_sizes;
        for (size_t i = 0; i < decomposition_sets_.size(); i++) {
            size_t real_size = 0;
            for (auto it = decomposition_sets_[i].begin(); it != decomposition_sets_[i].end(); it++)
                real_size += collapsed_struct_->GetMultiplicityOf(*it);
            real_sizes.push_back(real_size);
        }
        return real_sizes;
    }

    double AverageFillingInsideClass(const vector<size_t>& classes_filling) {
        vector<size_t> real_sizes = ComputeRealClassSizes();

        for(size_t i = 0; i < decomposition_sets_.size(); i++)
            if(real_sizes[i] > 1)
                cout << "Class: " << i << " size: " << decomposition_sets_[i].size() <<
                    ", num edges: " << classes_filling[i] << endl;

        ofstream out("decomposition_class_filling.txt", std::ios_base::app);
        double avg_filling = 0;
        size_t num_nt_classes = 0;
        for(size_t i = 0; i < decomposition_sets_.size(); i++)
            if(decomposition_sets_[i].size() > 1) {
                size_t class_size = decomposition_sets_[i].size();
                num_nt_classes++;
                double cur_filling = double(classes_filling[i]) /
                        double(class_size * (class_size - 1)) * 2;
                out << class_size << "\t" << cur_filling << "\t" << permutation_prt_->Size() << endl;
                avg_filling += cur_filling;
            }
        out.close();
        return avg_filling / double(num_nt_classes);
    }

    void ComputeStats() {
        double rel_num_classes = double(decomposition_sets_.size()) / double(permutation_prt_->Size());
        vector<size_t> classes_filling = ClassesFilling();
        double relative_class_filling = RelativeFilling(classes_filling);
        double avg_filling_inside_classes = AverageFillingInsideClass(classes_filling);
        TRACE("Decomposition statistics");
        TRACE("Relative number of classes: " << rel_num_classes);
        TRACE("% edges inside classes: " << relative_class_filling);
        TRACE("% edges outside classes: " << 1.0 - relative_class_filling);
        TRACE("Avg filling inside classes: " << avg_filling_inside_classes);

        ofstream out("decomposition_stats.txt", std::ios_base::app);
        out << rel_num_classes << "\t" << relative_class_filling << "\t" << 1 - relative_class_filling <<
                "\t" << avg_filling_inside_classes << endl;
    }
 */