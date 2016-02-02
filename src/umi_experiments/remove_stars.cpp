#include <logger/log_writers.hpp>
#include "../graph_utils/graph_io.hpp"

//using std::cout;
//using std::endl;
//using std::string;
using bformat = boost::format;

bool readArgs(int argc, char **argv, string& input_file, size_t& size_to_print) {
    if (argc < 2 || argc > 3) return false;
    input_file = argv[1];
    size_to_print = argc == 3 ? atoi(argv[2]) : 0;
    return true;
}

void create_console_logger() {
    using namespace logging;
    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

void mark(const SparseGraphPtr& graph, size_t v, vector<bool>& g, vector<size_t>& cc) {
    if (g[v]) return;
    g[v] = true;
    cc.push_back(v);
    for (auto u : graph->VertexEdges(v)) {
        mark(graph, u, g, cc);
    }
}

bool is_star(const SparseGraphPtr& graph, size_t v) {
    bool result = true;
    for (auto u : graph->VertexEdges(v)) {
        result &= graph->Weight()[u] < graph->Weight()[v] && graph->Degree(u) == 1;
    }
    return result;
}

bool is_doublet(const SparseGraphPtr& graph, size_t v) {
    return graph->Degree(v) == 1 && graph->Degree(*graph->VertexEdges(v).begin()) == 1;
}

int main(int argc, char **argv) {
    create_console_logger();
    string input_file = "";
    size_t size_to_print = 0;
    if (!readArgs(argc, argv, input_file, size_to_print)) {
        cout << "Usage: <input file> [<size of unclassified components to print>]" << endl;
        return 0;
    }
    const SparseGraphPtr& graph = GraphReader(input_file).CreateGraph();
    INFO("Searching and removing stars");
    vector<bool> g(graph->N(), false);
    size_t single = 0;
    size_t stars = 0;
    size_t stars_sum = 0;
    size_t doublets = 0;
    size_t left = 0;
    size_t left_components = 0;
    size_t gone = 0;
    for (size_t i = 0; i < graph->N(); i++) {
        if (g[i]) continue;
        vector<size_t> cc;
        mark(graph, i, g, cc);
        gone += cc.size();
        size_t size = cc.size();
        if (graph->Degree(i) == 0) {
            single++;
            continue;
        }
        if (graph->Degree(i) > 1) {
            if (is_star(graph, i)) {
                stars++;
                stars_sum += size;
                if (size < 2) {
                    ERROR(bformat("star of %d vertex") % size);
                }
                continue;
            }
        }
        {
            size_t v = *graph->VertexEdges(i).begin();
            if (graph->Weight()[i] > graph->Weight()[v]) {
                v = i;
            }
            if (is_star(graph, v)) {
                stars ++;
                stars_sum += size;
                continue;
            }
            if (is_doublet(graph, i)) {
                doublets ++;
                continue;
            }
        }
        left += size;
        left_components ++;
        if (size < size_to_print || size_to_print == 0) continue;
        stringstream ss;
        ss << endl << "graph" << endl;
        for (auto u : cc) {
            ss << graph->Weight()[u] << " ";
        }
        ss << endl;
        for (size_t u = 0; u < cc.size(); u ++) {
            for (size_t v = 0; v < u; v ++) {
                if (graph->HasEdge(cc[v], cc[u])) {
                    ss << "(" << v << ", " << u << ")  ";
                }
            }
        }
        INFO(ss.str());
    }
    INFO(bformat("Found %d singletons, %d stars, %d doublets, %d left, average star size %f, average other component size %f") % single % stars % doublets % left
         % (static_cast<double>(stars_sum) / static_cast<double>(stars)) % (static_cast<double>(left) / static_cast<double>(left_components)));
    assert(single + doublets * 2 + stars_sum + left == graph->N());
}
