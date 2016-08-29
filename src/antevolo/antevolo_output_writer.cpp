#include "antevolo_output_writer.hpp"

namespace antevolo {
    void AntEvoloOutputWriter::OutputSHMForTrees() const {
        std::ofstream out(output_params_.tree_shms);
        INFO("Tree SHMs were written to " << output_params_.tree_shms);
        out.close();
    }

    void AntEvoloOutputWriter::OutputTreeStats() const {
        std::ofstream out(output_params_.tree_details);
        out << "Tree_id\tNum_vertices\tSHM_depth\tRoot_depth\tNum_unique_SHMs\tNum_added_SHMs\tNum_synonymous_SHMs\t"
                       "Num_CDR3_SHMs\tNum_synonymous_CDR3_SHMs" << std::endl;
        size_t tree_index = 1;
        for(auto it = annotated_storage_.cbegin(); it != annotated_storage_.cend(); it++) {
            out << tree_index << "\t" << it->Tree().NumEdges() + 1 << "\t" << it->SHMDepth() << "\t" <<
                    it->RootDepth() << "\t" << it->NumUniqueSHms() << "\t" << it->NumAddedSHMs() <<
                    "\t" << it->NumSynonymousSHMs() << "\t" << it->NumCDR3SHMs() << "\t" << it->NumSynonymousCDR3SHMs() << std::endl;
            tree_index++;
//            if(it->SHMDepth() == 0) {
//                INFO("Tree with zero SHM depth:");
//                for(auto edge_it = it->Tree().cbegin(); edge_it != it->Tree().cend(); edge_it++) {
//                    std::cout << *edge_it << std::endl;
//                    std::cout << clone_set_[edge_it->src_clone_num].Read().seq << std::endl;
//                    std::cout << clone_set_[edge_it->dst_clone_num].Read().seq << std::endl;
//                }
//            }
        }
        INFO("Tree statistics were written to " << output_params_.tree_details);
        out.close();
    }

    void AntEvoloOutputWriter::WriteEdge(const EvolutionaryEdge& edge, std::ofstream& out) const { //no endl
        out << edge.src_clone->Read().id << "\t" << edge.dst_clone->Read().id << "\t"
            << edge.src_clone->Read().name << "\t" << edge.dst_clone->Read().name << "\t"
            << edge.src_clone->VSHMs().size() + edge.src_clone->JSHMs().size() << "\t"
            << edge.dst_clone->VSHMs().size() + edge.dst_clone->JSHMs().size()<< "\t"
            << edge.edge_type << "\t" <<  edge.num_intersected_shms << "\t" << edge.num_added_shms
            << "\t" << edge.cdr3_distance << "\t" << edge.weight << "\t"
            << edge.src_clone->Productive() << "\t" << edge.dst_clone->Productive() << "\t" << edge.IsSynonymous();
    }
    void AntEvoloOutputWriter::WriteTreeInFile(const EvolutionaryTree& tree) const {
        std::string output_fname = tree.GetTreeOutputFname();
        std::ofstream out(output_fname);
        out << "Src_id\tDst_id\tSrc_clone\tDst_clone\tNum_Src_SHMs\tNum_Dst_SHMs\tEdge_type\t";
        out << "Num_shared_SHMs\tNum_added_SHMs\tCDR3_dist\tWeight\t";
        out << "Src_productive\tDst_productive\tSynonymous\t";
        out << "Src_CDR3\tDst_CDR3\n";
        for(auto it = tree.cbegin(); it != tree.cend(); it++) {
            WriteEdge(*it, out);
            std::string src_CDR3_string;
            std::string dst_CDR3_string;
            auto const& src_CDR3 = it->src_clone->CDR3();
            auto const& dst_CDR3 = it->dst_clone->CDR3();
            size_t src_CDR3_length = seqan::length(src_CDR3);
            size_t dst_CDR3_length = seqan::length(dst_CDR3);
            for (size_t pos = 0; pos < src_CDR3_length; ++pos) {
                src_CDR3_string.push_back(src_CDR3[pos]);
            }
            for (size_t pos = 0; pos < dst_CDR3_length; ++pos) {
                dst_CDR3_string.push_back(dst_CDR3[pos]);
            }
            out << "\t" << src_CDR3_string << "\t" << dst_CDR3_string << std::endl;
        }
        out.close();
    }

    void AntEvoloOutputWriter::WriteTreeVerticesInFile(const EvolutionaryTree& tree,
                                               const annotation_utils::CDRAnnotatedCloneSet& clone_set) {
        std::string output_fname = tree.GetVerticesOutputFname();
        std::ofstream out(output_fname);
        out << "Clone_id\tClone_name\tProductive\tAA_seq\tOFR\tLeft_CDR3_anchor_AA\tRight_CDR3_anchor_AA\n";
        for (auto it = tree.c_vertex_begin(); it != tree.c_vertex_end(); it++) {
            auto const& clone = clone_set[*it];
            size_t ORF = clone.ORF();
            size_t start_pos = clone.GetRangeByRegion(
                    annotation_utils::StructuralRegion::CDR3).start_pos;
            size_t end_pos = clone.GetRangeByRegion(
                    annotation_utils::StructuralRegion::CDR3).end_pos;
            std::string clone_AA_string;
            auto const& clone_AA_seq = clone.AA();
            size_t AA_length = seqan::length(clone.AA());
            for (size_t i = 0; i < AA_length; ++i) {
                clone_AA_string.push_back(clone_AA_seq[i]);
            }

            //assert((static_cast<int>(start_pos) - static_cast<int>(ORF)) % 3 == 0);
            //assert((static_cast<int>(end_pos) + 1 - static_cast<int>(ORF)) % 3 == 0);
            char left_CDR3_anchor_AA = clone_AA_string[(start_pos - ORF) / 3 - 1];
            char right_CDR3_anchor_AA = clone_AA_string[(end_pos + 1 - ORF) / 3];

            out << clone.Read().id << "\t" << clone.Read().name << "\t" << clone.Productive() << "\t"
                << clone_AA_string << "\t" << ORF << "\t"
                << left_CDR3_anchor_AA << "\t" << right_CDR3_anchor_AA << "\n";

        }
        out.close();
    }
}