#include "antevolo_output_writer.hpp"
#include "shm_counting/tree_based_shm_convertor.hpp"

namespace antevolo {
    void AntEvoloOutputWriter::OutputSHMForTrees() const {
        for(auto it = annotated_storage_.cbegin(); it != annotated_storage_.cend(); it++) {
            std::string shm_fname = path::append_path(output_params_.tree_shm_dir, it->Tree().GetTreeOutputFname(""));
            std::ofstream out(shm_fname);
            WriteTreeSHMs(*it, out);
            out.close();
        }
        INFO("Tree SHMs were written to " << output_params_.tree_shm_dir);
    }

    void AntEvoloOutputWriter::WriteTreeSHMs(const AnnotatedEvolutionaryTree &tree, std::ofstream &out) const {
        auto shm_map = tree.SHMMap();
        auto root_id = tree.Tree().GetRoot();
        auto clone_set = tree.Tree().GetCloneSet();
        out << "@CDR1:" << clone_set[root_id].CDR1Range().start_pos << "," << clone_set[root_id].CDR1Range().end_pos << std::endl;
        out << "@CDR2:" << clone_set[root_id].CDR2Range().start_pos << "," << clone_set[root_id].CDR2Range().end_pos << std::endl;
        out << "@CDR3:" << clone_set[root_id].CDR3Range().start_pos << "," << clone_set[root_id].CDR3Range().end_pos << std::endl;
        out << "@VDJ_length:" << clone_set[root_id].Read().length() << std::endl;
        out << "Nucl_position\tAA_position\tSrc_nucl\tDst_nucl\tSrc_aa\tDst_aa\t"
                "Src_triplet\tDst_triplet\tMultiplicity\tRegion" << std::endl;
        std::map<size_t, std::vector<TreeSHM>> root_pos_shm_map;
        for(auto it = shm_map.c_shm_clone_begin(); it != shm_map.c_shm_clone_end(); it++) {
            auto shm = it->first;
//            auto shm_clone_id = it->second[0].second;
            auto root_n_pos = TreeSHMComparator::GetTreeSHMPosition(clone_set[root_id], shm);
            if(root_pos_shm_map.find(root_n_pos) == root_pos_shm_map.end()) {
                root_pos_shm_map[root_n_pos] = std::vector<TreeSHM>();
            }
            root_pos_shm_map[root_n_pos].push_back(shm);
        }
        for(auto it = root_pos_shm_map.begin(); it != root_pos_shm_map.end(); it++) {
            auto root_n_pos = it->first;
            auto shms = it->second;
            for(auto s = shms.begin(); s != shms.end(); s++) {
                auto shm = *s;
                out << root_n_pos << "\t" << clone_set[root_id].GetAminoAcidPosByNucleotidePos(root_n_pos) << "\t" <<
                    shm.src_nucl << "\t" << shm.dst_nucl << "\t" << shm.src_aa << "\t" << shm.dst_aa << "\t" <<
                    shm.src_triplet << "\t" << shm.dst_triplet << "\t" <<
                    it->second.size() << "\t" << shm.region << std::endl;
            }
        }
    }

    void AntEvoloOutputWriter::OutputTreeStats() const {
        std::ofstream out(output_params_.tree_details);
        out << "Tree_id\tNum_vertices\tTree_depth\tRoot_depth\tNum_unique_SHMs\tNum_synonymous_SHMs\t"
                    "Num_CDR_SHMs\tMax_multiplicity\t"
                    "FR1_len\tCDR1_len\tFR2_len\tCDR2_len\tFR3_len\tCDR3_len\tFR4_len\t"
                    "NumInFR1\tNumInCDR1\tNumInFR2\tNumInCDR2\t"
                    "NumInFR3\tNumInCDR3\tNumInFR4\t"
                    "NumSynInFR1\tNumSynInCDR1\tNumSynInFR2\tNumSynInCDR2\t"
                    "NumSynInFR3\tNumSynInCDR3\tNumSynInFR4" << std::endl;
        size_t tree_index = 1;
        for(auto it = annotated_storage_.cbegin(); it != annotated_storage_.cend(); it++) {
            out << tree_index << "\t" << it->Tree().NumVertices() << "\t" << it->TreeDepth() << "\t" <<
                it->RootDepth() << "\t" << it->NumUniqueSHMs() << "\t" <<
                it->NumSynonymousSHMs() << "\t" << it->NumSHMsInCDRs() << "\t" << it->MaxSHMMultiplicity() << "\t" <<
                //
                it->GetRegionLength(annotation_utils::StructuralRegion::FR1) << "\t" <<
                it->GetRegionLength(annotation_utils::StructuralRegion::CDR1) << "\t" <<
                it->GetRegionLength(annotation_utils::StructuralRegion::FR2) << "\t" <<
                it->GetRegionLength(annotation_utils::StructuralRegion::CDR2) << "\t" <<
                it->GetRegionLength(annotation_utils::StructuralRegion::FR3) << "\t" <<
                it->GetRegionLength(annotation_utils::StructuralRegion::CDR3) << "\t" <<
                it->GetRegionLength(annotation_utils::StructuralRegion::FR4) << "\t" <<
                //
                it->NumSHMsInRegion(annotation_utils::StructuralRegion::FR1) << "\t" <<
                it->NumSHMsInRegion(annotation_utils::StructuralRegion::CDR1) << "\t" <<
                it->NumSHMsInRegion(annotation_utils::StructuralRegion::FR2) << "\t" <<
                it->NumSHMsInRegion(annotation_utils::StructuralRegion::CDR2) << "\t" <<
                it->NumSHMsInRegion(annotation_utils::StructuralRegion::FR3) << "\t" <<
                it->NumSHMsInRegion(annotation_utils::StructuralRegion::CDR3) << "\t" <<
                it->NumSHMsInRegion(annotation_utils::StructuralRegion::FR4) << "\t" <<
                //
                it->NumSynonymousSHMsInRegion(annotation_utils::StructuralRegion::FR1) << "\t" <<
                it->NumSynonymousSHMsInRegion(annotation_utils::StructuralRegion::CDR1) << "\t" <<
                it->NumSynonymousSHMsInRegion(annotation_utils::StructuralRegion::FR2) << "\t" <<
                it->NumSynonymousSHMsInRegion(annotation_utils::StructuralRegion::CDR2) << "\t" <<
                it->NumSynonymousSHMsInRegion(annotation_utils::StructuralRegion::FR3) << "\t" <<
                it->NumSynonymousSHMsInRegion(annotation_utils::StructuralRegion::CDR3) << "\t" <<
                it->NumSynonymousSHMsInRegion(annotation_utils::StructuralRegion::FR4) << std::endl;
            tree_index++;
        }
        INFO("Clonal tree statistics were written to " << output_params_.tree_details);
        out.close();
    }

    void AntEvoloOutputWriter::WriteEdge(const EvolutionaryEdgePtr& edge, std::ofstream& out) const { //no endl
        out << edge->SrcClone()->Read().id << "\t" << edge->DstClone()->Read().id << "\t"
            << edge->SrcClone()->Read().name << "\t" << edge->DstClone()->Read().name << "\t"
            << edge->SrcClone()->VSHMs().size() + edge->SrcClone()->JSHMs().size() << "\t"
            << edge->DstClone()->VSHMs().size() + edge->DstClone()->JSHMs().size()<< "\t"
            << edge->TypeString() << "\t" <<  edge->NumSharedShms() << "\t" << edge->NumAddedShms()
            << "\t" << edge->CDR3Distance() << "\t" << edge->Length() << "\t"
            << edge->SrcClone()->Productive() << "\t" << edge->DstClone()->Productive() << "\t" << edge->IsSynonymous();
    }
    void AntEvoloOutputWriter::WriteTreeInFile(std::string output_dir, const EvolutionaryTree& tree) const {
        std::string output_fname = tree.GetTreeOutputFname(output_dir);
        std::ofstream out(output_fname);
        out << "Src_id\tDst_id\tSrc_clone\tDst_clone\tNum_Src_SHMs\tNum_Dst_SHMs\tEdge_type\t";
        out << "Num_shared_SHMs\tNum_added_SHMs\tCDR3_dist\tWeight\t";
        out << "Src_productive\tDst_productive\tSynonymous\t";
        out << "Src_CDR3\tDst_CDR3\n";
        for(auto it = tree.cbegin(); it != tree.cend(); it++) {
            WriteEdge(*it, out);
            std::string src_CDR3_string;
            std::string dst_CDR3_string;
            auto const& src_CDR3 = (*it)->SrcClone()->CDR3();
            auto const& dst_CDR3 = (*it)->DstClone()->CDR3();
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

    void AntEvoloOutputWriter::WriteTreeVerticesInFile(std::string output_dir, const EvolutionaryTree& tree) const {
        const auto& clone_set = tree.GetCloneSet();
        std::string output_fname = tree.GetTreeOutputFname(output_dir);
        std::ofstream out(output_fname);
        out << "Clone_id\tClone_name\tProductive\tAA_seq\tOFR\tLeft_CDR3_anchor_AA\tRight_CDR3_anchor_AA\tSize\tFake\n";
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

            //VERIFY((static_cast<int>(start_pos) - static_cast<int>(ORF)) % 3 == 0);
            //VERIFY((static_cast<int>(end_pos) + 1 - static_cast<int>(ORF)) % 3 == 0);
            char left_CDR3_anchor_AA = clone_AA_string[(start_pos - ORF) / 3 - 1];
            char right_CDR3_anchor_AA = clone_AA_string[(end_pos + 1 - ORF) / 3];

            out << clone.Read().id << "\t" << clone.Read().name << "\t" << clone.Productive()
                << "\t" << clone_AA_string << "\t" << ORF << "\t"
                << left_CDR3_anchor_AA << "\t" << right_CDR3_anchor_AA << "\t" << clone.Size()
                << "\t" << tree.GetCloneSet().IsFake(*it) << "\n";

        }
        out.close();
    }

    void AntEvoloOutputWriter::OutputCleanedSequences(CloneSetWithFakesPtr clone_set_ptr) const {
        const auto& clone_set = *clone_set_ptr;
        std::string output_fname = path::append_path(output_params_.output_dir, "cleaned_sequences.fa");
        std::ofstream out(output_fname);
        for (size_t i = 0; i < clone_set.size(); ++i) {
            const auto& clone = clone_set[i];
            out << ">" << clone.Read().name << "\n";
            out << clone.Read().seq << "\n";
        }
        out.close();
    }

    void AntEvoloOutputWriter::WriteRcmFromStorageInFile(std::string output_dir,
                                                         const EvolutionaryTreeStorage& storage) {

        if (storage.size() == 0) {
            return;
        }
        const auto& clone_set = storage[0].GetCloneSet().GetOriginalCloneSet();
        std::vector<bool> written_down(clone_set.size());
        std::string output_fname = path::append_path(output_dir, "clonal_lineage_decomposition.rcm");
        std::ofstream out(output_fname);
        size_t current_cluster;
        for (current_cluster = 0; current_cluster < storage.size(); ++current_cluster) {
            const auto& tree = storage[current_cluster];
            for (auto it = tree.c_vertex_begin(); it != tree.c_vertex_end(); ++it) {
                if (tree.GetCloneSet().IsFake(*it)) {
                    continue;
                }
                out << tree.GetCloneSet()[*it].Read().name << "\t" << current_cluster << "\n";
                written_down[*it] = true;
            }
        }
        for (size_t i = 0; i < clone_set.size(); ++i) {
            if (!written_down[i]) {
                out << clone_set[i].Read().name << "\t" << current_cluster << "\n";
                ++current_cluster;
                written_down[i] = true;
            }
        }
        out.close();
    }
}