#include <verify.hpp>

#include "diversity_analyser.hpp"
#include "../graph_utils/graph_splitter.hpp"

#include <../graph_utils/graph_io.hpp>
#include <../graph_utils/graph_splitter.hpp>

namespace cdr_labeler {
    size_t DiversityAnalyser::ComputeD50(const std::vector<SparseGraphPtr> connected_components) const {
        std::vector<size_t> sizes;
        for(auto it = connected_components.begin(); it != connected_components.end(); it++)
            sizes.push_back((*it)->N());
        std::sort(sizes.begin(), sizes.end());
        size_t sum = 0;
        for(auto it = sizes.begin(); it != sizes.end(); it++) {
            sum += *it;
            if(sum >= cdr_graph_->N() / 2)
                return *it;
        }
        return size_t(-1);
    }

    void DiversityAnalyser::InitializeGraph(std::string compressed_cdr3_fasta) {
        std::string graph_fname = "cdr3.graph";
        std::string run_graph_constructor = "./build/release/bin/ig_swgraph_construct";
        std::string command_line = run_graph_constructor + " -i " + compressed_cdr3_fasta +
                " -o " + graph_fname + " --tau 3 > out";
        int err_code = system(command_line.c_str());
        VERIFY_MSG(err_code == 0, "Graph constructor finished abnormally, error code: " << err_code);
        cdr_graph_ = GraphReader(graph_fname).CreateGraph();
        ConnectedComponentGraphSplitter graph_splitter(cdr_graph_);
        auto connected_components = graph_splitter.Split();
        INFO(connected_components.size() << " connected components of CDR3 Hamming graph were created");
        size_t num_isolated = 0;
        size_t size_largest = 0;
        for(auto it = connected_components.begin(); it != connected_components.end(); it++) {
            if((*it)->N() == 1)
                num_isolated++;
            size_largest = std::max<size_t>(size_largest, (*it)->N());
        }
        INFO("# isolated CDR3s: " << num_isolated);
        INFO("Largest component size: " << size_largest);
        INFO("D50: " << ComputeD50(connected_components));
    }
}