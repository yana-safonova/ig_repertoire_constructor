#include "crs_matrix.hpp"

class GraphReader {
    std::string graph_filename;

public:
    GraphReader(std::string graph_filename) {
        this->graph_filename = graph_filename;
    }

    CRS_HammingGraph_Ptr CreateGraph();
};

class GraphWriter {
    std::string graph_filename;

public:
    GraphWriter(std::string graph_filename) {
        this->graph_filename = graph_filename;
    }

    void PrintGraph(CRS_HammingGraph_Ptr graph_ptr);
};