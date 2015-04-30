#pragma once
#include "../../include_me.hpp"

/*
    struct HGEdge characterizes edge of Hamming graph
    contains fields:
        i & j - incident vertices,
        dist - hamming distance between these vertices
        shm - number of SHM patterns these vertices
 */
struct HGEdge {
    size_t i;
    size_t j;
    size_t dist;
    size_t shm;

    HGEdge(size_t i, size_t j, size_t dist, size_t shm = 0):
            i(i),
            j(j),
            dist(dist),
            shm(shm) { }
};

/*
    class CRS_Matrix for storage of sparse symmetric matrix in compressed rows format
 */
class CRS_Matrix {
    size_t N_;
    size_t NZ_;

    // rows and columns
    vector<size_t> row_index_;
    vector<size_t> col_;

    // weights
    vector<size_t> dist_;

    void InitializeRowIndex() {
        for(size_t i = 0; i < N_+ 1; i++)
            row_index_.push_back(0);
    }

    // valid only for consecutive filling of structure
    void AddEdge(size_t v1, size_t v2, size_t dist) {
        assert(v1 < v2);
        dist_.push_back(dist);
        col_.push_back(v2);
        row_index_[v1]++;
    }

    void ComputeRowIndex() {
        size_t num_next_elem = 0;
        for(size_t i = 0; i < N_ + 1; i++) {
            size_t tmp = row_index_[i];
            row_index_[i] = num_next_elem;
            num_next_elem += tmp;
        }
    }

    void ComputeN_NZ(const vector<HGEdge> &edges) {
        N_ = 0;
        for(auto e = edges.begin(); e != edges.end(); e++)
            N_ = max<size_t>(N_, max<size_t>(e->i, e->j));
        N_++;
        NZ_ = edges.size();
    }

    // edges are consecutive
    void Initialize(const vector<HGEdge> &edges) {
        ComputeN_NZ(edges);
        InitializeRowIndex();
        for(auto it = edges.begin(); it != edges.end(); it++)
            AddEdge(it->i, it->j, it->dist);
        ComputeRowIndex();
    }

    // NOTE: input matrix is transposed
    void Initialize(const CRS_Matrix &trans_matrix) {
        N_ = trans_matrix.N();
        NZ_ = trans_matrix.NZ();
        vector<vector<size_t> > trans_nz_col;
        vector<vector<size_t> > trans_nz_val;
        for(size_t i = 0; i < trans_matrix.N(); i++) {
            trans_nz_col.push_back(vector<size_t>());
            trans_nz_val.push_back(vector<size_t>());
        }
        for(size_t i = 0; i < trans_matrix.N(); i++)
            for(size_t j = trans_matrix.RowIndex()[i]; j < trans_matrix.RowIndex()[i + 1]; j++) {
                size_t row = i;
                size_t col = trans_matrix.Col()[j];
                size_t value = trans_matrix.Dist()[j];
                trans_nz_col[col].push_back(row);
                trans_nz_val[col].push_back(value);
            }
        row_index_.push_back(0);
        for(size_t i = 0; i < trans_nz_col.size(); i++) {
            for(size_t j = 0; j < trans_nz_col[i].size(); j++) {
                col_.push_back(trans_nz_col[i][j]);
                dist_.push_back(trans_nz_val[i][j]);
            }
            row_index_.push_back(row_index_[row_index_.size() - 1] + trans_nz_col[i].size());
        }
    }

public:
    CRS_Matrix(const CRS_Matrix &trans_matrix) :
            N_(0), NZ_(0) {
        Initialize(trans_matrix);
    }

    CRS_Matrix(const vector<HGEdge> &edges) :
            N_(0), NZ_(0) {
        Initialize(edges);
    }

    const vector<size_t>& Dist() const { return dist_; }

    const vector<size_t>& RowIndex() const  { return row_index_; }

    const vector<size_t>& Col() const { return col_; }

    size_t N() const { return N_; }

    size_t NZ() const { return NZ_; }

    shared_ptr<CRS_Matrix> Transpose() {
        return shared_ptr<CRS_Matrix>(new CRS_Matrix(*this));
    }
};

ostream& operator<<(ostream &out, const CRS_Matrix &matrix) {
    out << "N: " << matrix.N() << ", NZ: " << matrix.NZ() << endl;
    out << "Columns: ";
    for(auto it = matrix.Col().begin(); it != matrix.Col().end(); it++)
        out << *it << " ";
    out << endl << "Values: ";
    for(auto it = matrix.Dist().begin(); it != matrix.Dist().end(); it++)
        out << *it << " ";
    out << endl << "Row indices: ";
    for(auto it = matrix.RowIndex().begin(); it != matrix.RowIndex().end(); it++)
        out << *it << " ";
    return out;
}

typedef shared_ptr<CRS_Matrix> CRS_Matrix_Ptr;

/*
    class Hamming graph in CRS format
    contains upper and lower triangle representation of graph matrix
 */
class CRS_HammingGraph {
    CRS_Matrix_Ptr direct_matrix_;
    CRS_Matrix_Ptr trans_matrix_;

public:
    CRS_HammingGraph(const vector<HGEdge> &edges) :
            direct_matrix_(new CRS_Matrix(edges)) {
        trans_matrix_ = direct_matrix_->Transpose();
    }

    size_t N() const { return direct_matrix_->N(); }

    size_t NZ() const { return direct_matrix_->NZ(); }

    size_t Degree(size_t i) const  {
        assert(i < N());
        return direct_matrix_->RowIndex()[i + 1] - direct_matrix_->RowIndex()[i] +
                trans_matrix_->RowIndex()[i + 1] - trans_matrix_->RowIndex()[i];
    }

    const vector<size_t>& RowIndex() const { return direct_matrix_->RowIndex(); }

    const vector<size_t>& RowIndexT() const { return trans_matrix_->RowIndex(); }

    const vector<size_t>& Col() const { return direct_matrix_->Col(); }

    const vector<size_t>& ColT() const { return trans_matrix_->Col(); }

    const vector<size_t>& Dist() const { return direct_matrix_->Dist(); }

    const vector<size_t>& DistT() const { return trans_matrix_->Dist(); }

    const CRS_Matrix_Ptr DirectMatrix() const { return direct_matrix_; }

    const CRS_Matrix_Ptr TransposedMatrix() const { return trans_matrix_; }
};

ostream& operator<<(ostream &out, const CRS_HammingGraph &graph) {
    out << "Direct matrix" << endl;
    out << *(graph.DirectMatrix()) << endl;
    out << "-------------" << endl;
    out << "Transposed matrix" << endl;
    out << *(graph.TransposedMatrix());
    return out;
}

typedef shared_ptr<CRS_HammingGraph> CRS_HammingGraph_Ptr;

class HG_CollapsedStructs {
    CRS_HammingGraph_Ptr hamming_graph_ptr_;

    size_t num_main_vertices_;
    size_t num_all_vertices_;
    vector<size_t> old_vertices_list_; // old_vertices_list[i] shows old number of new vertex i
    vector<size_t> new_vertices_list_; // new_vertices_list[i] shows new number of old vertex i
    vector<size_t> main_vertices_map_; // main_vertices_map[i] shows index of main vertex (old indices)
    vector<size_t> multiplicity_;

    void InitializeMainVerticesMap() {
        vector<size_t> main_vertices_tree;
        for(size_t i = 0; i < hamming_graph_ptr_->N(); i++) {
            multiplicity_.push_back(1);
            main_vertices_tree.push_back(i);
        }

        for(size_t i = 0; i < hamming_graph_ptr_->N(); i++)
            for(size_t j = hamming_graph_ptr_->RowIndex()[i]; j < hamming_graph_ptr_->RowIndex()[i + 1];j++)
                if(hamming_graph_ptr_->Dist()[j] == 0)
                    main_vertices_tree[hamming_graph_ptr_->Col()[j]] = i;

        for(size_t i = 0; i < main_vertices_tree.size(); i++) {
            size_t cur_vertex = i;
            while(cur_vertex != main_vertices_tree[cur_vertex])
                cur_vertex = main_vertices_tree[cur_vertex];
            main_vertices_map_.push_back(cur_vertex);
            multiplicity_[cur_vertex]++;
        }
    }

    void InitializeOldVerticesList() {
        for(size_t i = 0; i < main_vertices_map_.size(); i++)
            if(main_vertices_map_[i] == i)
                old_vertices_list_.push_back(i);
    }

    void InitializeNewVerticesList() {
        for(size_t i = 0; i < hamming_graph_ptr_->N(); i++)
            new_vertices_list_.push_back(size_t(-1));
        for(size_t i = 0; i < old_vertices_list_.size(); i++)
            new_vertices_list_[old_vertices_list_[i]] = i;
    }

    void Initialize() {
        InitializeMainVerticesMap();
        InitializeOldVerticesList();
        InitializeNewVerticesList();
    }

public:
    HG_CollapsedStructs(CRS_HammingGraph_Ptr hamming_graph_ptr) :
            hamming_graph_ptr_(hamming_graph_ptr) {
        Initialize();
    }

    const vector<size_t>& MainVerticesMap() const { return main_vertices_map_; }

    const vector<size_t>& OldVerticesList() const { return old_vertices_list_; }

    const vector<size_t>& NewVerticesList() const { return new_vertices_list_; }

    size_t NumCollapsedVertices() const { return old_vertices_list_.size(); }

    bool OldIndexValid(size_t index) {
        assert(index < new_vertices_list_.size());
        return new_vertices_list_[index] != size_t(-1);
    }

    bool NewIndexValid(size_t index) {
        return index != size_t(-1);
    }

    size_t GetMultiplicityOf(size_t index) {
        assert(index < multiplicity_.size());
        return multiplicity_[index];
    }

    bool VertexIsMain(size_t index) {
        assert(index < main_vertices_map_.size());
        return index == main_vertices_map_[index];
    }

    size_t NewIndexOfOldVertex(size_t old_vertex) {
        assert(old_vertex < main_vertices_map_.size());
        return new_vertices_list_[main_vertices_map_[old_vertex]];
    }

    size_t NumberNewVertices() { return old_vertices_list_.size(); }

    size_t MultiplicityOfNewVertex(size_t new_index) {
        assert(new_index < old_vertices_list_.size());
        return multiplicity_[old_vertices_list_[new_index]];
    }

    size_t MultiplicityOfOldVertex(size_t old_vertex) {
        assert(old_vertex < multiplicity_.size());
        return multiplicity_[old_vertex];
    }
};

typedef shared_ptr<HG_CollapsedStructs> HG_CollapsedStructs_Ptr;