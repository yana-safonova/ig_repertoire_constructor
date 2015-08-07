#include "graph_io.hpp"
#include "../../ig_tools/utils/string_tools.hpp"

/*
 * class WeightedGraphReader
 *      takes as an input file stream starting from the first line of the adjacency list
 *      creates list of weighted edges
 *      returns Hamming graph
 */
class WeightedGraphReader {
    vector<GraphEdge> graph_edges;

    void UpdateGraphEdges(size_t cur_vertex, const vector<string> &line_splits) {
        assert(line_splits.size() % 2 == 0);
        size_t num_edges = line_splits.size() / 2;
        for(size_t i = 0; i < num_edges; i++) {
            size_t dst_vertex = string_to_number<size_t>(line_splits[i * 2]);
            size_t edge_weight = string_to_number<size_t>(line_splits[i * 2 + 1]);
            graph_edges.push_back(GraphEdge(cur_vertex, dst_vertex, edge_weight));
        }
    }

public:
    WeightedGraphReader() { }

    CRS_HammingGraph_Ptr ReadGraph(std::ifstream &graph_stream) {
        size_t cur_vertex = 0;
        while(!graph_stream.eof()) {
            std::string tmp_line;
            getline(graph_stream, tmp_line);
            vector<string> splits = split(tmp_line, ' ');
            UpdateGraphEdges(cur_vertex, splits);
            cur_vertex++;
        }
        return CRS_HammingGraph_Ptr(new CRS_HammingGraph(cur_vertex, graph_edges));
    }
};

/*
 * class UnweightedGraphReader
 *      takes as an input file stream starting from the first line of the adjacency list
 *      creates list of edges of weight 1
 *      returns Hamming graph
 */
class UnweightedGraphReader {
    vector<GraphEdge> graph_edges;

    void UpdateGraphEdges(size_t cur_vertex, const vector<string> &line_splits) {
        size_t num_edges = line_splits.size();
        for(size_t i = 0; i < num_edges; i++) {
            size_t dst_vertex = string_to_number<size_t>(line_splits[i]);
            size_t edge_weight = 1;
            graph_edges.push_back(GraphEdge(cur_vertex, dst_vertex, edge_weight));
        }
    }

public:
    UnweightedGraphReader() { }

    CRS_HammingGraph_Ptr ReadGraph(std::ifstream &graph_stream) {
        size_t cur_vertex = 0;
        while(!graph_stream.eof()) {
            std::string tmp_line;
            getline(graph_stream, tmp_line);
            vector<string> splits = split(tmp_line, ' ');
            UpdateGraphEdges(cur_vertex, splits);
            cur_vertex++;
        }
        return CRS_HammingGraph_Ptr(new CRS_HammingGraph(cur_vertex, graph_edges));
    }
};

/*
 *
 */
class VersatileGraphReader {
    bool GraphIsUnweighted(const vector<string> &header_splits) {
        return header_splits.size() == 2;
    }

    bool GraphInWeighted(const vector<string> &header_splits) {
        if(header_splits.size() != 3)
            return false;
        return header_splits[2] == "001";
    }

public:
    CRS_HammingGraph_Ptr ReadGraph(std::ifstream &graph_stream) {
        string header_line;
        getline(graph_stream, header_line);
        vector<string> splits = split(header_line, ' ');
        if(GraphIsUnweighted(splits))
            return UnweightedGraphReader().ReadGraph(graph_stream);
        if(GraphInWeighted(splits))
            return WeightedGraphReader().ReadGraph(graph_stream);
        return WeightedGraphReader().ReadGraph(graph_stream);
    }
};

/*
 *
 */
CRS_HammingGraph_Ptr GraphReader::CreateGraph() {
    std::ifstream graph_stream(this->graph_filename);
    if(graph_stream.fail()) {
        std::cout << "File " + this->graph_filename + " with graph was not found" << std::endl;
    }
    CRS_HammingGraph_Ptr graph_ptr = VersatileGraphReader().ReadGraph(graph_stream);
    INFO("Graph contains " << graph_ptr->N() << " vertices and " << graph_ptr->NZ() << " edges");
    return graph_ptr;
}