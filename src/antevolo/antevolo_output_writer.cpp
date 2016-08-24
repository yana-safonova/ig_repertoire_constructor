#include "antevolo_output_writer.hpp"

namespace antevolo {
    void AntEvoloOutputWriter::OutputSHMForTrees() const {
        std::ofstream out(output_params_.tree_shms);
        INFO("Tree SHMs were written to " << output_params_.tree_shms);
        out.close();
    }

    void AntEvoloOutputWriter::OutputTreeStats() const {
        std::ofstream out(output_params_.tree_details);
        out << "Tree_id\tNum_vertices\tSHM_depth\tRoot_depth\tNum_unique_SHMs\tNum_added_SHMs\tNum_synonymous_SHMs" << std::endl;
        size_t tree_index = 1;
        for(auto it = annotated_storage_.cbegin(); it != annotated_storage_.cend(); it++) {
            out << tree_index << "\t" << it->Tree().NumEdges() + 1 << "\t" << it->SHMDepth() << "\t" <<
                    it->RootDepth() << "\t" << it->NumUniqueSHms() << "\t" << it->NumAddedSHMs() <<
                    "\t" << it->NumSynonymousSHMs() << std::endl;
            tree_index++;
        }
        INFO("Tree statistics were written to " << output_params_.tree_details);
        out.close();
    }
}