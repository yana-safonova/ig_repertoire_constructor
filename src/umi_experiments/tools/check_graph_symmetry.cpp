#include <logger/log_writers.hpp>
#include "../graph_utils/sparse_graph.hpp"
#include "../graph_utils/graph_io.hpp"
#include "utils.hpp"

int main(int, char **argv) {
    create_console_logger();
    string input_file = argv[1];
    const SparseGraphPtr graph = GraphReader(input_file).CreateGraph();
    vector<vector<size_t> > g(graph->N());
    int e = 0;
    for (size_t v = 0; v < graph->N(); v ++) {
        g[v] = vector<size_t>();
        for (size_t u : graph->VertexEdges(v)) {
            g[v].push_back(u);
            e ++;
        }
    }
    INFO(boost::format("Read graph with v = %d, e = %d") % graph->N() % graph->NZ());
    for (size_t v = 0; v < g.size(); v ++) {
        for (auto u: g[v]) {
            bool found = false;
            for (auto w: g[u]) {
                if (w == v) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                ERROR(boost::format("There's an edge %d->%d, but no %d->%d") % v % u % u % v);
                return 0;
            }
        }
    }
}