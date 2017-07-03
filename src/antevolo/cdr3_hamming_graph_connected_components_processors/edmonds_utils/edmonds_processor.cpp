#include "edmonds_processor.hpp"
#include "../../antevolo_launch.hpp"

namespace antevolo {

    std::vector<WeightedEdge<int>> EdmondsProcessor::process_edge_list(
            const std::vector<WeightedEdge<int>>& input_edges)
    {
        if (input_edges.size() == 0) {
            return std::vector<WeightedEdge<int>>();
        }
//        std::cout << "here it begins" << std::endl;
        size_t n = 0;
        boost::unordered_map<size_t, size_t> vertex_to_index;
        std::vector<size_t> index_to_vertex;
        for (auto e : input_edges) {
            if (vertex_to_index.find(e.src_) == vertex_to_index.end()) {
                vertex_to_index[e.src_] = n;
                index_to_vertex.push_back(e.src_);
                ++n;
            }
            if (vertex_to_index.find(e.dst_) == vertex_to_index.end()) {
                vertex_to_index[e.dst_] = n;
                index_to_vertex.push_back(e.dst_);
                ++n;
            }
        }

        Graph G(n);
//        for (auto v : index_to_vertex) {
//            std::cout << v << " ";
//        } std::cout << std::endl;
//
//        for (auto v : vertex_to_index) {
//            std::cout << v.first << ": " << v.second << " | ";
//        } std::cout << std::endl;
//        std::cout << "here it calculates maps" << std::endl;
        std::vector<Vertex> the_vertices;
        BOOST_FOREACH (Vertex v, vertices(G))
        {
            the_vertices.push_back(v);
        }

        for (auto e : input_edges) {
            size_t src_index = vertex_to_index[e.src_];
            size_t dst_index = vertex_to_index[e.dst_];
            add_edge(the_vertices[src_index], the_vertices[dst_index], e.weight_, G);
        }

//        std::cout << "here it adds edges" << std::endl;

        boost::property_map<Graph, boost::edge_weight_t>::type weights =
                get(boost::edge_weight_t(), G);

        boost::property_map<Graph, boost::vertex_index_t>::type vertex_indices =
                get(boost::vertex_index_t(), G);



        std::vector<Vertex> roots = {the_vertices[0]};
        std::vector<Edge> branching;
        edmonds_optimum_branching<false, true, false>(G,
                                                     vertex_indices,
                                                     weights,
                                                     static_cast<Vertex *>(0),
                                                     static_cast<Vertex *>(0),
                                                     std::back_inserter(branching));
//        std::cout << "here it calculates branching" << std::endl;
        std::vector<WeightedEdge<int>> res;
//        INFO(n);
//        for (auto v : index_to_vertex) {
//            std::cout << v << " ";
//        } std::cout << std::endl;
        BOOST_FOREACH (Edge e, branching)
        {
//            INFO(boost::source(e, G) << " " << boost::target(e, G) << " " << boost::source(e, G));
            size_t src = index_to_vertex[boost::source(e, G)];
            size_t dst = index_to_vertex[boost::target(e, G)];
            int weight = get(weights, e);
            res.push_back(WeightedEdge<int>(src, dst, weight));
        }
//        std::cout << "here it ends" << std::endl;
        return res;
    }

}