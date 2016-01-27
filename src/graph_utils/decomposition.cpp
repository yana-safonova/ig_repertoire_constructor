#include "decomposition.hpp"
#include "../ig_tools/utils/string_tools.hpp"

Decomposition::Decomposition(string decomposition_filename) {
    ifstream in(decomposition_filename);
    assert(in.good());
    vector <size_t> classes_list = ReadClassIdsFromIfstream(in);
    TRACE("Decomposition of size " << classes_list.size() << " was extracted from " << decomposition_filename);
    num_vertices_ = classes_list.size();
    InitializeVertexClasses();
    for (size_t i = 0; i < classes_list.size(); i++)
        SetClass(i, classes_list[i]);
}

void Decomposition::InitializeVertexClasses() {
    num_classes_ = 0;
    for(size_t i = 0; i < num_vertices_; i++)
        vertex_class_.push_back(size_t(-1));
}

vector<size_t> Decomposition::ReadClassIdsFromIfstream(ifstream &in) {
    vector<size_t> classes_list;
    while(!in.eof()) {
        string tmp;
        getline(in, tmp);
        if(tmp == "")
            break;
        classes_list.push_back(string_to_number<size_t>(tmp));
    }
    return classes_list;
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

void Decomposition::AddDecomposition(shared_ptr<Decomposition> decomposition) {
    size_t last_class_id = LastClassId();
    for(size_t i = 0; i < decomposition->Size(); i++) {
        set<size_t> cur_class = decomposition->GetClass(i);
        for(auto it = cur_class.begin(); it != cur_class.end(); it++)
            SetClass(*it, last_class_id + i);
    }
};

void Decomposition::SaveTo(string output_fname) {
    ofstream out(output_fname.c_str());
    for(auto it = vertex_class_.begin(); it != vertex_class_.end(); it++)
        out << *it << endl;
    out.close();
}

size_t Decomposition::MaxClassSize() {
    size_t max_class = 0;
    for(auto it = decomposition_classes_.begin(); it != decomposition_classes_.end(); it++)
        max_class = max<size_t>(max_class, it->size());
    return max_class;
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