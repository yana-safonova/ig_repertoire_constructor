#include "clonal_graph_writer.hpp"

namespace antevolo {
    void ClonalGraphWriter::operator()(std::string output_fname) {
        std::ofstream out(output_fname);
        out << "digraph G {" << std::endl;
        for(auto e = cgraph_.cbegin(); e != cgraph_.cend(); e++)  {
            auto src = e->first;
            auto dst_vertices = e->second;
            for(auto it = dst_vertices.begin(); it != dst_vertices.end(); it++) {
                out << src << " -> " << *it << std::endl;
            }
        }
        out << "}" << std::endl;
        out.close();
    }
}