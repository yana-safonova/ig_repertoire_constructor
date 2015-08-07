#include "include_me.hpp"

#include "graph_utils/crs_matrix.hpp"
#include "graph_utils/graph_io.hpp"

#include "metis_permutation_constructor/metis_permutation_constructor.hpp"

int main(int argc, char* argv[]) {
    if(argc != 2) {
        std::cout << "dense_sgraph_finder graph.GRAPH" << std::endl;
        return 1;
    }
    std::string graph_fname = std::string(argv[1]);
    GraphReader graph_reader(graph_fname);
    CRS_HammingGraph_Ptr graph_ptr = graph_reader.CreateGraph();

    //MetisPermutationConstructor permutation_constructor();

}