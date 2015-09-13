#include "crs_matrix.hpp"

// valid only for consecutive filling of structure
void CrsMatrix::AddEdge(size_t v1, size_t v2, size_t dist) {
    assert(v1 < v2);
    dist_.push_back(dist);
    col_.push_back(v2);
    row_index_[v1]++;
}

void CrsMatrix::ComputeRowIndex() {
    size_t num_next_elem = 0;
    for(size_t i = 0; i < N_ + 1; i++) {
        size_t tmp = row_index_[i];
        row_index_[i] = num_next_elem;
        num_next_elem += tmp;
    }
}

void CrsMatrix::Initialize(const vector<GraphEdge> &edges) {
    ComputeN_NZ(edges);
    InitializeRowIndex();
    for(auto it = edges.begin(); it != edges.end(); it++)
        AddEdge(it->i, it->j, it->dist);
    ComputeRowIndex();
}

// NOTE: input matrix is transposed
void CrsMatrix::Initialize(const CrsMatrix &trans_matrix) {
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

ostream& operator<<(ostream &out, const CrsMatrix &matrix) {
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