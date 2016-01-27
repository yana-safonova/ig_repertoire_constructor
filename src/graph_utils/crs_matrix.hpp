#pragma once
#include "include_me.hpp"

/*
    struct Edge characterizes edge of graph
    contains fields:
        i & j - incident vertices,
        dist - edge weight
 */
struct GraphEdge {
    size_t i;
    size_t j;
    size_t dist;

    GraphEdge(size_t i, size_t j, size_t dist):
            i(i),
            j(j),
            dist(dist) { }
};

/*
    class CRS_Matrix for storage of sparse symmetric matrix in compressed rows format
 */
class CrsMatrix {
    size_t N_;
    size_t NZ_;

    // rows and columns
    vector<size_t> row_index_;
    vector<size_t> col_;

    // weights
    vector<size_t> dist_;

    // isolated vertices
    vector<size_t> isolated_vertices_;

    void InitializeRowIndex() {
        for(size_t i = 0; i < N_+ 1; i++)
            row_index_.push_back(0);
    }

    // valid only for consecutive filling of structure
    void AddEdge(size_t v1, size_t v2, size_t dist);

    void ComputeRowIndex();

    void ComputeN_NZ(const vector<GraphEdge> &edges) {
        NZ_ = edges.size();
    }

    // edges are consecutive
    void Initialize(const vector<GraphEdge> &edges);

    // NOTE: input matrix is transposed
    void Initialize(const CrsMatrix &trans_matrix);

public:
    CrsMatrix(const CrsMatrix &trans_matrix) :
            N_(0), NZ_(0) {
        Initialize(trans_matrix);
    }

    CrsMatrix(size_t N, const vector<GraphEdge> &edges) :
            N_(N), NZ_(0) {
        Initialize(edges);
    }

    const vector<size_t>& Dist() const { return dist_; }

    const vector<size_t>& RowIndex() const  { return row_index_; }

    const vector<size_t>& Col() const { return col_; }

    size_t N() const { return N_; }

    size_t NZ() const { return NZ_; }

    std::shared_ptr<CrsMatrix> Transpose() {
        return std::shared_ptr<CrsMatrix>(new CrsMatrix(*this));
    }
};

ostream& operator<<(ostream &out, const CrsMatrix &matrix);

typedef std::shared_ptr<CrsMatrix> CrsMatrixPtr;