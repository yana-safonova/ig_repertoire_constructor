#pragma once

#include "graph_collapsed_structure.hpp"

class Decomposition {
    // input params
    size_t num_vertices_;

    // decomposition fields
    vector<size_t> vertex_class_;
    vector<set<size_t> > decomposition_classes_;

    // number of all classes: real and removed
    size_t num_classes_;

    void InitializeVertexClasses();

    vector<size_t> ReadClassIdsFromIfstream(ifstream &in);

    void AddNewClass();

    bool ClassIsValid(size_t class_id) { return class_id != size_t(-1); }

    void RemoveVertex(size_t vertex);

public:
    Decomposition(size_t num_vertices) :
            num_vertices_(num_vertices) {
        InitializeVertexClasses();
    }

    Decomposition(string decomposition_filename);

    void SetClass(size_t vertex, size_t class_id);

    void AddDecomposition(shared_ptr<Decomposition> decomposition);

    const set<size_t>& GetClass(size_t index) const {
        assert(index < Size());
        return decomposition_classes_[index];
    }

    const set<size_t>& LastClass() const { return GetClass(Size() - 1); }

    size_t ClassSize(size_t index) const { return GetClass(index).size(); }

    size_t LastClassSize() const { return LastClass().size(); }

    size_t LastClassId() const { return Size() - 1; }

    size_t VertexNumber() const { return vertex_class_.size(); }

    size_t GetVertexClass(size_t vertex) const {
        assert(vertex < VertexNumber());
        return vertex_class_[vertex];
    }

    bool VertexClassIsInitialized(size_t vertex) {
        return ClassIsValid(GetVertexClass(vertex));
    }

    size_t Size() const { return decomposition_classes_.size(); }

    void SaveTo(string output_fname);

    bool LastClassContains(size_t vertex) const {
        return LastClass().find(vertex) != LastClass().end();
    }

    size_t MaxClassSize();

    bool IsTrivial() { return Size() == 1; }

    size_t NextClassId() { return num_classes_; }

private:
    DECL_LOGGER("Decomposition");
};

typedef set<size_t> DecompositionClass;

ostream& operator<<(ostream &out, const Decomposition &hg_decomposition);

typedef shared_ptr<Decomposition> DecompositionPtr;
