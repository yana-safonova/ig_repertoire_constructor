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
