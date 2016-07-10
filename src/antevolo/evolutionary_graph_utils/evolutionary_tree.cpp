#include "evolutionary_tree.hpp"

namespace antevolo {
    void EvolutionaryTree::WriteInFile(std::string output_fname) {
        std::ofstream out(output_fname);
        out << "Src_num\tDst_num\tSrc_clone\tDst_clone\tNum_Src_V_SHMs\tNum_Dst_V_SHMs\tEdge_type\tNum_shared_SHMs\tNum_added_SHMs\tCDR3_dist\tWeight\n";
        for(auto it = edges_.begin(); it != edges_.end(); it++)
            out << it->second.src_clone_num << "\t" << it->second.dst_clone_num << "\t" << it->second.src_clone->Read().name << "\t" << it->second.dst_clone->Read().name << "\t" <<
            it->second.src_clone->VSHMs().size() << "\t" << it->second.dst_clone->VSHMs().size() << "\t" << it->second.edge_type << "\t" <<
            it->second.num_intersected_v_shms << "\t" << it->second.num_added_v_shms << "\t" << it->second.cdr3_distance << "\t" <<
            it->second.weight << std::endl;
        out.close();
    }
}