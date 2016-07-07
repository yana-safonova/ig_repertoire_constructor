#include "evolutionary_tree.hpp"

namespace antevolo {
    void EvolutionaryTree::WriteInFile(std::string output_fname) {
        std::ofstream out(output_fname);
        out << "Src_clone\tDst_clone\tNum_Src_V_SHMs\tNum_Dst_V_SHMs\tEdge_type\tNum_shared_SHMs\tNum_added_SHMs\tCDR3_dist\tWeight\n";
        for(auto it = edges_.begin(); it != edges_.end(); it++)
            out << it->src_clone->Read().name << "\t" << it->dst_clone->Read().name << "\t" <<
            it->src_clone->VSHMs().size() << "\t" << it->dst_clone->VSHMs().size() << "\t" << it->edge_type << "\t" <<
            it->num_intersected_v_shms << "\t" << it->num_added_v_shms << "\t" << it->cdr3_distance << "\t" <<
            it->weight << std::endl;
        out.close();
    }
}