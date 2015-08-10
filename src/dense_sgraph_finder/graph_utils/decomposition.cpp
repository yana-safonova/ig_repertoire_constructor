#include "decomposition.hpp"

void Decomposition::InitializeVertexClasses() {
    num_classes_ = 0;
    for(size_t i = 0; i < num_vertices_; i++)
        vertex_class_.push_back(size_t(-1));
}

void Decomposition::AddNewClass() {
    set<size_t> new_class;
    decomposition_classes_.push_back(new_class);
    num_classes_++;
}

void Decomposition::RemoveVertex(size_t vertex) {
    size_t vertex_class = vertex_class_[vertex];
    if(!ClassIsValid(vertex_class))
        return;
    decomposition_classes_[vertex_class].erase(vertex);
}

void Decomposition::SetClass(size_t vertex, size_t class_id) {
    assert(vertex < num_vertices_);
    RemoveVertex(vertex);
    if(decomposition_classes_.size() == 0 or class_id > decomposition_classes_.size() - 1)
        for(size_t i = decomposition_classes_.size(); i <= class_id; i++)
            AddNewClass();
    decomposition_classes_[class_id].insert(vertex);
    vertex_class_[vertex] = class_id;
}

void Decomposition::SaveTo(string output_fname) {
    ofstream out(output_fname.c_str());
    for(auto it = vertex_class_.begin(); it != vertex_class_.end(); it++)
        out << *it << endl;
    out.close();
    INFO("Decomposition was written to " << output_fname);
}

size_t Decomposition::MaxClassSize() {
    size_t max_class = 0;
    for(auto it = decomposition_classes_.begin(); it != decomposition_classes_.end(); it++)
        max_class = max<size_t>(max_class, it->size());
    return max_class;
}

size_t Decomposition::RealSizeOfClass(size_t class_id, GraphCollapsedStructurePtr collapsed_struct_ptr) {
    assert(class_id < Size());
    size_t class_size = 0;
    for(auto it = decomposition_classes_[class_id].begin();
        it != decomposition_classes_[class_id].end(); it++) {
        size_t old_index = collapsed_struct_ptr->OldVerticesList()[*it];
        class_size += collapsed_struct_ptr->GetMultiplicityOf(old_index);
    }
    return class_size;
}

ostream& operator<<(ostream &out, const Decomposition &hg_decomposition) {
    out << "Decomposition contains " << hg_decomposition.Size() << " classes" << endl;
    for(size_t i = 0; i < hg_decomposition.Size(); i++) {
        out << "Class " << i << ", size: " << hg_decomposition.ClassSize(i) << ". ";
        for(auto it = hg_decomposition.GetClass(i).begin(); it != hg_decomposition.GetClass(i).end(); it++)
            out << *it << " ";
        out << endl;
    }
    return out;
}