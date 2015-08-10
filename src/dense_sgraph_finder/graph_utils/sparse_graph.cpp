#include "sparse_graph.hpp"

ostream& operator<<(ostream &out, const SparseGraph &graph) {
    out << "Direct matrix" << endl;
    out << *(graph.DirectMatrix()) << endl;
    out << "-------------" << endl;
    out << "Transposed matrix" << endl;
    out << *(graph.TransposedMatrix());
    return out;
}