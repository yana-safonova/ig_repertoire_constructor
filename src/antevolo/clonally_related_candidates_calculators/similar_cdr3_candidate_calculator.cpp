#include <path_helper.hpp>
#include "similar_cdr3_candidate_calculator.hpp"
#include "../../graph_utils/graph_io.hpp"
#include "../../graph_utils/graph_splitter.hpp"
#include <convert.hpp>
#include <annotation_utils/shm_comparator.hpp>
#include <evolutionary_graph_utils/evolutionary_graph_constructor.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <map>

namespace antevolo {
    void SimilarCDR3CandidateCalculator::Clear() {
        unique_cdr3s_.clear();
        unique_cdr3s_map_.clear();
    }

    void SimilarCDR3CandidateCalculator::CreateUniqueCDR3Map(
            core::DecompositionClass decomposition_class) {
        for(auto it = decomposition_class.begin(); it != decomposition_class.end(); it++) {
            if(clone_set_[*it].RegionIsEmpty(annotation_utils::StructuralRegion::CDR3))
                continue;
            auto cdr3 = core::dna5String_to_string(clone_set_[*it].CDR3());
            if(unique_cdr3s_map_.find(cdr3) == unique_cdr3s_map_.end())
                unique_cdr3s_map_[cdr3] = std::vector<size_t>();
            unique_cdr3s_map_[cdr3].push_back(*it);
        }
        for(auto it = unique_cdr3s_map_.begin(); it != unique_cdr3s_map_.end(); it++)
            unique_cdr3s_.push_back(it->first);
    }

    std::string SimilarCDR3CandidateCalculator::GetFastaFname(core::DecompositionClass decomposition_class) {
        std::stringstream ss;
        size_t key = *decomposition_class.begin();
        ss << "CDR3_sequences_key_" << key << ".fasta";
        return path::append_path(output_params_.cdr_graph_dir, ss.str());
    }

    std::string SimilarCDR3CandidateCalculator::GetGraphFname(core::DecompositionClass decomposition_class) {
        std::stringstream ss;
        size_t key = *decomposition_class.begin();
        ss << "CDR3_sequences_key_" << key << ".graph";
        return path::append_path(output_params_.cdr_graph_dir, ss.str());
    }

    std::string SimilarCDR3CandidateCalculator::WriteUniqueCDR3InFasta(core::DecompositionClass decomposition_class) {
        std::string output_fname = GetFastaFname(decomposition_class);
        std::ofstream out(output_fname);
        for(auto it = unique_cdr3s_.begin(); it != unique_cdr3s_.end(); it++)
            out << ">" << *it << std::endl << *it << std::endl;
        return output_fname;
    }

    // return connected components of Hamming graph on CDR3s
    std::vector<SparseGraphPtr> SimilarCDR3CandidateCalculator::ComputeCDR3HammingGraphs(std::string cdr_fasta,
                                                                                         std::string graph_fname) {
        std::string run_graph_constructor = "./build/release/bin/ig_swgraph_construct";
        std::stringstream ss;
        ss << run_graph_constructor << " -i " << cdr_fasta <<
                " -o " << graph_fname << " --tau " << num_mismatches_ << " -k 10 > " << output_params_.trash_output;
        int err_code = system(ss.str().c_str());
        VERIFY_MSG(err_code == 0, "Graph constructor finished abnormally, error code: " << err_code);
        auto sparse_cdr_graph_ = GraphReader(graph_fname).CreateGraph();
        TRACE("Hamming graph contains " << sparse_cdr_graph_->N() << " edges and " << sparse_cdr_graph_->NZ() << " edges");
        //if(sparse_cdr_graph_->NZ() != 0)
        //    std::cout << *sparse_cdr_graph_ << std::endl;
        auto connected_components = ConnectedComponentGraphSplitter(sparse_cdr_graph_).Split();
        graph_component_ = sparse_cdr_graph_->GetGraphComponentMap();
        return connected_components;
    }

    /*
    ClonallyRelatedCandidates SimilarCDR3CandidateCalculator::ComputeCandidatesForGraph(SparseGraphPtr hg_component,
                                                                                        size_t component_id) {
        ClonallyRelatedCandidates candidates(clone_set_);
        // adding edges between identical CDR3s
        for(size_t i = 0; i < hg_component->N(); i++) {
            size_t old_index = graph_component_.GetOldVertexByNewVertex(component_id, i);
            auto clones_sharing_cdr3 = unique_cdr3s_map_[unique_cdr3s_[old_index]];
            for(size_t it1 = 0; it1 < clones_sharing_cdr3.size(); it1++)
                for(size_t it2 = it1 + 1; it2 < clones_sharing_cdr3.size(); it2++)
                    candidates.AddCandidatePair(clones_sharing_cdr3[it1], clones_sharing_cdr3[it2]);
        }
        // adding edges between similar CDR3s
        for(size_t i = 0; i < hg_component->N(); i++)
            for(size_t j = hg_component->RowIndex()[i]; j < hg_component->RowIndex()[i + 1]; j++) {
                size_t old_index1 = graph_component_.GetOldVertexByNewVertex(component_id, i);
                size_t old_index2 = graph_component_.GetOldVertexByNewVertex(component_id, hg_component->Col()[j]);
                auto indices_1 = unique_cdr3s_map_[unique_cdr3s_[old_index1]];
                auto indices_2 = unique_cdr3s_map_[unique_cdr3s_[old_index2]];
                for(auto it1 = indices_1.begin(); it1!= indices_1.end(); it1++)
                    for(auto it2 = indices_2.begin(); it2!= indices_2.end(); it2++)
                        candidates.AddCandidatePair(*it1, *it2);
            }
        return candidates;
    }
     */
    /*
    EvolutionaryTree SimilarCDR3CandidateCalculator::ComputeTreeForComponent(SparseGraphPtr hg_component,
                                                                                        size_t component_id) {
        EvolutionaryTree tree;
        auto edge_constructor = GetEdgeConstructor();
        // adding edges between identical CDR3s
        for(size_t i = 0; i < hg_component->N(); i++) {
            size_t old_index = graph_component_.GetOldVertexByNewVertex(component_id, i);
            auto clones_sharing_cdr3 = unique_cdr3s_map_[unique_cdr3s_[old_index]];
            for(size_t it1 = 0; it1 < clones_sharing_cdr3.size(); it1++)
                for(size_t it2 = it1 + 1; it2 < clones_sharing_cdr3.size(); it2++) {
                    auto edge = edge_constructor->ConstructEdge(clone_set_[it1], clone_set_[it2], it1, it2);
                    tree.Add(clones_sharing_cdr3[it2], edge);
                }
        }
        // adding edges between similar CDR3s
        for(size_t i = 0; i < hg_component->N(); i++)
            for(size_t j = hg_component->RowIndex()[i]; j < hg_component->RowIndex()[i + 1]; j++) {
                size_t old_index1 = graph_component_.GetOldVertexByNewVertex(component_id, i);
                size_t old_index2 = graph_component_.GetOldVertexByNewVertex(component_id, hg_component->Col()[j]);
                auto indices_1 = unique_cdr3s_map_[unique_cdr3s_[old_index1]];
                auto indices_2 = unique_cdr3s_map_[unique_cdr3s_[old_index2]];
                for(auto it1 = indices_1.begin(); it1!= indices_1.end(); it1++)
                    for(auto it2 = indices_2.begin(); it2!= indices_2.end(); it2++) {
                        auto edge = edge_constructor->ConstructEdge(clone_set_[*it1], clone_set_[*it2], *it1, *it2);
                        tree.Add(*it2, edge);
                    }
            }
        return tree;
    }
    */

    void SimilarCDR3CandidateCalculator::AddComponentToTheTree(SparseGraphPtr hg_component,
                                                                             size_t component_id, EvolutionaryTree& tree) {
        auto edge_constructor = GetEdgeConstructor();

        std::vector<size_t> vertices_nums;
        for(size_t i = 0; i < hg_component->N(); i++) {
            vertices_nums.push_back(graph_component_.GetOldVertexByNewVertex(component_id, i));
        }
        typedef boost::associative_property_map<std::map<size_t, size_t>> AP_map;
        AP_map rank;
        AP_map parent;
        boost::disjoint_sets<AP_map, AP_map> ds_on_undirected_edges(rank, parent);
        for (size_t i : vertices_nums) {
            ds_on_undirected_edges.make_set(i);
        }



        //adding undirected edges first
        for(size_t i = 0; i < hg_component->N(); i++) {
            size_t old_index = graph_component_.GetOldVertexByNewVertex(component_id, i);
            auto clones_sharing_cdr3 = unique_cdr3s_map_[unique_cdr3s_[old_index]];
            for(size_t it1 = 0; it1 < clones_sharing_cdr3.size(); it1++)
                for(size_t it2 = it1 + 1; it2 < clones_sharing_cdr3.size(); it2++) {
                    size_t src_num = clones_sharing_cdr3[it1];
                    size_t dst_num = clones_sharing_cdr3[it2];
                    // if clonea are not in the same connected component and
                    // the edge is undirected
                    if (ds_on_undirected_edges.find_set(src_num) !=
                            ds_on_undirected_edges.find_set(dst_num) &&
                            annotation_utils::SHMComparator::SHMsAreEqual(
                                    clone_set_[src_num].VSHMs(),
                                    clone_set_[dst_num].VSHMs())) {
                        tree.AddUndirectedPair(src_num, dst_num);
                        ds_on_undirected_edges.union_set(src_num, dst_num);
                    }
                }
        }
        for(size_t i = 0; i < hg_component->N(); i++)
            for(size_t j = hg_component->RowIndex()[i]; j < hg_component->RowIndex()[i + 1]; j++) {
                size_t old_index1 = graph_component_.GetOldVertexByNewVertex(component_id, i);
                size_t old_index2 = graph_component_.GetOldVertexByNewVertex(component_id, hg_component->Col()[j]);
                auto indices_1 = unique_cdr3s_map_[unique_cdr3s_[old_index1]];
                auto indices_2 = unique_cdr3s_map_[unique_cdr3s_[old_index2]];
                for(auto it1 = indices_1.begin(); it1!= indices_1.end(); it1++)
                    for(auto it2 = indices_2.begin(); it2!= indices_2.end(); it2++) {
                        if (ds_on_undirected_edges.find_set(*it1) !=
                                    ds_on_undirected_edges.find_set(*it2) &&
                                annotation_utils::SHMComparator::SHMsAreEqual(
                                        clone_set_[*it1].VSHMs(),
                                        clone_set_[*it2].VSHMs())) {
                            tree.AddUndirectedPair(*it1, *it2);
                            ds_on_undirected_edges.union_set(*it1, *it2);
                        }
                    }
            }

        // adding directed edges between identical CDR3s
        for(size_t i = 0; i < hg_component->N(); i++) {
            size_t old_index = graph_component_.GetOldVertexByNewVertex(component_id, i);
            auto clones_sharing_cdr3 = unique_cdr3s_map_[unique_cdr3s_[old_index]];
            for(size_t it1 = 0; it1 < clones_sharing_cdr3.size(); it1++)
                for(size_t it2 = it1 + 1; it2 < clones_sharing_cdr3.size(); it2++) {
                    size_t src_num = clones_sharing_cdr3[it1];
                    size_t dst_num = clones_sharing_cdr3[it2];
                    auto edge = edge_constructor->ConstructEdge(
                            clone_set_[src_num],
                            clone_set_[dst_num],
                            src_num,
                            dst_num);
                    tree.AddDirected(dst_num, edge);
                    std::vector<std::pair<size_t, size_t>> edge_vector;
                    tree.Prepare_subtree(edge_vector, dst_num);
                    for (auto p : edge_vector) {
                        auto edge = edge_constructor->ConstructEdge(
                                clone_set_[p.first],
                                clone_set_[p.second],
                                p.first,
                                p.second);
                        tree.AddUndirected(p.second, edge);
                    }
                }
        }
        // adding directed edges between similar CDR3s
        for(size_t i = 0; i < hg_component->N(); i++)
            for(size_t j = hg_component->RowIndex()[i]; j < hg_component->RowIndex()[i + 1]; j++) {
                size_t old_index1 = graph_component_.GetOldVertexByNewVertex(component_id, i);
                size_t old_index2 = graph_component_.GetOldVertexByNewVertex(component_id, hg_component->Col()[j]);
                auto indices_1 = unique_cdr3s_map_[unique_cdr3s_[old_index1]];
                auto indices_2 = unique_cdr3s_map_[unique_cdr3s_[old_index2]];
                for(auto it1 = indices_1.begin(); it1!= indices_1.end(); it1++)
                    for(auto it2 = indices_2.begin(); it2!= indices_2.end(); it2++) {
                        auto edge = edge_constructor->ConstructEdge(
                                clone_set_[*it1],
                                clone_set_[*it2],
                                *it1,
                                *it2);
                        tree.AddDirected(*it2, edge);
                        std::vector<std::pair<size_t, size_t>> edge_vector;
                        tree.Prepare_subtree(edge_vector, *it2);
                        for (auto p : edge_vector) {
                            auto edge = edge_constructor->ConstructEdge(
                                    clone_set_[p.first],
                                    clone_set_[p.second],
                                    p.first,
                                    p.second);
                            tree.AddUndirected(p.second, edge);
                        }
                    }
            }
        auto undirected_graph = tree.Get_undirected_graph();
        for (auto p : undirected_graph) {
            std::vector<std::pair<size_t, size_t>> edge_vector;
            tree.Prepare_subtree(edge_vector, p.first);
            for (auto p : edge_vector) {
                auto edge = edge_constructor->ConstructEdge(
                        clone_set_[p.first],
                        clone_set_[p.second],
                        p.first,
                        p.second);
                tree.AddUndirected(p.second, edge);
            }
        }



    }

    /*
    std::vector<ClonallyRelatedCandidates> SimilarCDR3CandidateCalculator::ComputeCandidates(
            core::DecompositionClass decomposition_class) {
        Clear();
        CreateUniqueCDR3Map(decomposition_class);
        std::string cdrs_fasta = WriteUniqueCDR3InFasta(decomposition_class);
        std::string graph_fname = GetGraphFname(decomposition_class);
        TRACE("--------------------------");
        TRACE("CDR3 fasta: "<< cdrs_fasta << ", CDR3 Hamming graph: " << graph_fname);
        auto connected_components = ComputeCDR3HammingGraphs(cdrs_fasta, graph_fname);
        TRACE("# connected components: " << connected_components.size());
        std::vector<ClonallyRelatedCandidates> vector_candidates;
        for(size_t i = 0; i < connected_components.size(); i++)
            vector_candidates.push_back(ComputeCandidatesForGraph(connected_components[i], i));
        return vector_candidates;
    }
     */


    EvolutionaryEdgeConstructor* SimilarCDR3CandidateCalculator::GetEdgeConstructor() {
        return new SimpleEvolutionaryEdgeConstructor(config_.edge_construction_params);
    }


    /*
    std::vector<std::pair<size_t, size_t>> SimilarCDR3CandidateCalculator::ComputeDirectedEdges(
            core::DecompositionClass decomposition_class) {
        Clear();
        CreateUniqueCDR3Map(decomposition_class);
        std::pair<size_t, std::string> top_cdr(0, "");
        std::pair<size_t, std::string> top2_cdr(0, "");
        for (auto it = unique_cdr3s_.begin(); it != unique_cdr3s_.end(); it++) {
            size_t cur_cdr_size = unique_cdr3s_map_[*it].size();
            if (cur_cdr_size > top2_cdr.first) {
                if (cur_cdr_size > top_cdr.first) {
                    top2_cdr = top_cdr;
                    top_cdr = std::make_pair(cur_cdr_size, *it);
                }
            }
        }
        std::vector<size_t> top2_cdr_set = unique_cdr3s_map_[top2_cdr.second];
        std::vector<std::pair<size_t, size_t>> top2_cdr_clones;
        for (auto clone_num : unique_cdr3s_map_[top2_cdr.second]) {
            top2_cdr_clones.push_back(std::make_pair(clone_set_[clone_num].VSHMs().size(), clone_num));
        }
        //sort by the VSHms number
        std::sort(top2_cdr_clones.begin(), top2_cdr_clones.end());
        std::vector<std::pair<size_t, size_t>> edges_to_add;
        for (auto it1 = top2_cdr_clones.rbegin(); it1 != top2_cdr_clones.rend(); it1++) {
            for (auto it2 = it1 + 1; it2 != top2_cdr_clones.rend(); it2++) {
                if (annotation_utils::SHMComparator::SHMs1AreNestedInSHMs2(clone_set_[it1->second].VSHMs(),
                                                                           clone_set_[it2->second].VSHMs())) {
                    edges_to_add.push_back(std::make_pair(it2->second, it1->second));
                    break;
                }
            }
        }
        return edges_to_add;
    }
     */
}