#include <verify.hpp>
#include "graph_io.hpp"
#include "../ig_tools/utils/string_tools.hpp"

vector<string> SplitGraphString(string str) {
    vector<string> final_splits;
    vector<string> splits = split(str, ' ');
    for(auto it = splits.begin(); it != splits.end(); it++) {
        vector<string> cur_splits = split(*it, '\t');
        for(auto it2 = cur_splits.begin(); it2 != cur_splits.end(); it2++)
            final_splits.push_back(*it2);
    }
    return final_splits;
}

/*
 * class EdgeWeightedGraphReader
 *      takes as an input file stream starting from the first line of the adjacency list
 *      creates list of weighted edges
 *      returns graph
 */
class EdgeWeightedGraphReader {
    vector<GraphEdge> graph_edges;

    void UpdateGraphEdges(size_t cur_vertex, const vector<string> &line_splits) {
        if(line_splits.size() % 2 != 0) {
            WARN("Line for vertex " << cur_vertex << " contains odd number of elements (" << line_splits.size() << "):");
            for(auto it = line_splits.begin(); it != line_splits.end(); it++)
                cout << *it << " ";
            cout << endl;
            assert(line_splits.size() % 2 == 0);
        }
        for(size_t i = 0; i < line_splits.size(); i += 2) {
            size_t dst_vertex = string_to_number<size_t>(line_splits[i]) - 1;
            size_t edge_weight = string_to_number<size_t>(line_splits[i + 1]);
            if(cur_vertex < dst_vertex)
                graph_edges.push_back(GraphEdge(cur_vertex, dst_vertex, edge_weight));
        }
    }

public:
    EdgeWeightedGraphReader() { }

    SparseGraphPtr ReadGraph(size_t num_vertices, ifstream &graph_stream) {
        size_t cur_vertex = 0;
        while(!graph_stream.eof()) {
            std::string tmp_line;
            getline(graph_stream, tmp_line);
            vector<string> splits = SplitGraphString(tmp_line);
            if (splits.size() == 0) {
                cur_vertex++;
                continue;
            }
            UpdateGraphEdges(cur_vertex, splits);
            cur_vertex++;
        }
        return SparseGraphPtr(new SparseGraph(num_vertices, graph_edges));
    }
};

/*
 * class EdgeVertexWeightedGraphReader
 *      takes as an input file stream starting from the first line of the adjacency list
 *      creates list of weighted edges and weighted vertices
 *      returns graph
 */
class EdgeVertexWeightedGraphReader {
    vector<GraphEdge> graph_edges;

    void UpdateGraphEdges(size_t cur_vertex, const vector<string> &line_splits) {
        if(line_splits.size() % 2 == 0) {
            WARN("Line for vertex " << cur_vertex << " contains even number of elements (" << line_splits.size() << "):");
            for(auto it = line_splits.begin(); it != line_splits.end(); it++)
                cout << *it << " ";
            cout << endl;
            assert(line_splits.size() % 2 != 0);
        }
        for(size_t i = 1; i < line_splits.size(); i += 2) {
            size_t dst_vertex = string_to_number<size_t>(line_splits[i]) - 1;
            size_t edge_weight = string_to_number<size_t>(line_splits[i + 1]);
            if(cur_vertex < dst_vertex)
                graph_edges.push_back(GraphEdge(cur_vertex, dst_vertex, edge_weight));
        }
    }

public:
    EdgeVertexWeightedGraphReader() { }

    SparseGraphPtr ReadGraph(size_t num_vertices, ifstream &graph_stream) {
        size_t cur_vertex = 0;
        auto weight = vector<size_t>(num_vertices, 1);
        while(!graph_stream.eof()) {
            std::string tmp_line;
            getline(graph_stream, tmp_line);
            vector<string> splits = SplitGraphString(tmp_line);
            if (splits.size() == 0)
                continue;
            weight[cur_vertex] = string_to_number<size_t>(splits[0]);
            UpdateGraphEdges(cur_vertex, splits);
            cur_vertex++;
        }
        return SparseGraphPtr(new SparseGraph(num_vertices, graph_edges, weight));
    }
};

/*
 * class UnweightedGraphReader
 *      takes as an input file stream starting from the first line of the adjacency list
 *      creates list of edges of weight 1
 *      returns graph
 */
class UnweightedGraphReader {
    vector<GraphEdge> graph_edges;

    void UpdateGraphEdges(size_t cur_vertex, const vector<string> &line_splits) {
        size_t num_edges = line_splits.size();
        for(size_t i = 0; i < num_edges; i++) {
            size_t dst_vertex = string_to_number<size_t>(line_splits[i]) - 1;
            size_t edge_weight = 1;
            if(cur_vertex < dst_vertex)
                graph_edges.push_back(GraphEdge(cur_vertex, dst_vertex, edge_weight));
        }
    }

public:
    UnweightedGraphReader() { }

    SparseGraphPtr ReadGraph(size_t num_vertices, std::ifstream &graph_stream) {
        size_t cur_vertex = 0;
        while(!graph_stream.eof()) {
            std::string tmp_line;
            getline(graph_stream, tmp_line);
            vector<string> splits = SplitGraphString(tmp_line);
            UpdateGraphEdges(cur_vertex, splits);
            cur_vertex++;
        }
        return SparseGraphPtr(new SparseGraph(num_vertices, graph_edges));
    }
};

/*
 *
 */
class VersatileGraphReader {
    bool GraphIsUnweighted(const vector<string> &header_splits) {
        return header_splits.size() == 2;
    }

    bool GraphIsEdgeWeighted(const vector <string> &header_splits) {
        if(header_splits.size() != 3 || header_splits[2].length() != 3)
            return false;
        return header_splits[2] == "001";
    }

    bool GraphIsEdgeVertexWeighted(const vector <string> &header_splits) {
        if(header_splits.size() != 3 || header_splits[2].length() != 3)
            return false;
        return header_splits[2] == "011";
    }

    size_t GetNumVertices(const vector<string> &header_splits) {
        assert(header_splits.size() > 0);
        return string_to_number<size_t>(header_splits[0]);
    }

public:
    SparseGraphPtr ReadGraph(std::ifstream &graph_stream) {
        string header_line;
        getline(graph_stream, header_line);
        vector<string> splits = SplitGraphString(header_line);
        size_t num_vertices = GetNumVertices(splits);
        if (GraphIsUnweighted(splits)) {
            TRACE("Unweighted graph reader was chosen");
            return UnweightedGraphReader().ReadGraph(num_vertices, graph_stream);
        }
        if (GraphIsEdgeWeighted(splits)) {
            TRACE("Edge weighted graph reader was chosen");
            return EdgeWeightedGraphReader().ReadGraph(num_vertices, graph_stream);
        }
        if (GraphIsEdgeVertexWeighted(splits)) {
            TRACE("Edge & vertex graph reader was chosen");
            return EdgeVertexWeightedGraphReader().ReadGraph(num_vertices, graph_stream);
        }
        VERIFY_MSG(false, "Unknown format of graph file!");
        return SparseGraphPtr(NULL);
    }

private:
    DECL_LOGGER("VersatileGraphReader");
};

/*
 *
 */
SparseGraphPtr GraphReader::CreateGraph() {
    TRACE("Trying to extract graph from " + graph_filename);
    std::ifstream graph_stream(graph_filename);
    if(!graph_stream.good()) {
        WARN("File " + this->graph_filename + " with graph was not found");
        return SparseGraphPtr(NULL);
    }
    SparseGraphPtr graph_ptr = VersatileGraphReader().ReadGraph(graph_stream);
    TRACE("Extracted graph contains " << graph_ptr->N() << " vertices & " << graph_ptr->NZ() << " edges");
    return graph_ptr;
}