#pragma once

#include "../graph_utils/sparse_graph.hpp"
#include "../graph_utils/graph_io.hpp"
#include "graph_decomposer/dense_subgraph_constructor.hpp"

namespace dense_subgraph_finder {
    class DenseSubgraphFinder {
        const dsf_config::run_params &run_params_;
        const dsf_config::dense_sgraph_finder_params &dsf_params_;
        const dsf_config::io_params &io_;
        const dsf_config::metis_io_params &metis_io_;

    public:
        DenseSubgraphFinder(const dsf_config::run_params &run_params,
                            const dsf_config::dense_sgraph_finder_params &dsf_params,
                            const dsf_config::io_params &io,
                            const dsf_config::metis_io_params &metis_io) :
                run_params_(run_params),
                dsf_params_(dsf_params),
                io_(io),
                metis_io_(metis_io) { }

        int Run();
    };
}