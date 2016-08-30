#include <verify.hpp>

#include "model_utils/shm_model.hpp"
#include <clonally_related_candidates_calculators/edmonds_tarjan_DMST_calculator.hpp>
#include "evolutionary_tree.hpp"

namespace antevolo {
    void EvolutionaryTree::AddEdge(size_t dst_id, EvolutionaryEdge edge) {
        VERIFY(dst_id == edge.dst_clone_num);
        edges_[dst_id] = edge;
        all_edge_vector_.push_back(edge);
        vertices_.insert(edge.dst_clone_num);
        vertices_.insert(edge.src_clone_num);
        if(outgoing_edges_.find(edge.src_clone_num) == outgoing_edges_.end()) {
            outgoing_edges_[edge.src_clone_num] = std::vector<EvolutionaryEdge>();
        }
        outgoing_edges_[edge.src_clone_num].push_back(edge);
    }

    void EvolutionaryTree::AddDirected(size_t clone_num, EvolutionaryEdge edge) {
        VERIFY(edge.IsDirected());
        //std::cout << "oppa: " << clone_num << " - " << edge.src_clone_num << " - " << edge.dst_clone_num << std::endl;
        if(edge.IsDirected()) {
            if (!Contains(clone_num)) {
                AddEdge(clone_num, edge);
                return;
            }
            const EvolutionaryEdge& parent_edge =  edges_[clone_num];
            if (parent_edge.num_added_shms > edge.num_added_shms) {
                //if clone_set_[*it2] is root or if the new edge is shorter
                AddEdge(clone_num, edge);
                return;
            }
            if (parent_edge.num_added_shms == edge.num_added_shms && parent_edge.cdr3_distance > edge.cdr3_distance) {
                //if clone_set_[*it2] is root or if the new edge is shorter
                AddEdge(clone_num, edge);
                return;
            }
        }
    }

    /*
    void EvolutionaryTree::SetUndirectedComponentParentEdge(size_t root_num, EvolutionaryEdge edge, ShmModel& model) {
        if(edge.IsDirected()) {
            if (undirected_components_edges_.find(root_num) == undirected_components_edges_.end()) {
                undirected_components_edges_[root_num] = edge;
                return;
            }
            const EvolutionaryEdge& parent_edge =  undirected_components_edges_[root_num];
            if (parent_edge.num_added_shms > edge.num_added_shms) {
                //if clone_set_[*it2] is root or if the new edge is shorter
                undirected_components_edges_[root_num] = edge;
                return;
            }
            if (parent_edge.num_added_shms == edge.num_added_shms && parent_edge.cdr3_distance > edge.cdr3_distance) {
                //if clone_set_[*it2] is root or if the new edge is shorter
                undirected_components_edges_[root_num] = edge;
                return;
            }
        }
    }
    */

    void EvolutionaryTree::AddUndirected(size_t clone_num, EvolutionaryEdge edge) {
        VERIFY(edge.IsUndirected());
        if (edge.IsUndirected()) {
            AddEdge(clone_num, edge);
        }
    }

    /*
    void EvolutionaryTree::AddUndirectedPair(size_t src_num, size_t dst_num) {
        if (undirected_graph_.find(src_num) == undirected_graph_.end()) {
            undirected_graph_[src_num] = std::set<size_t>();
        }
        if (undirected_graph_.find(dst_num) == undirected_graph_.end()) {
            undirected_graph_[dst_num] = std::set<size_t>();
        }
        if (undirected_graph_[src_num].find(dst_num) == undirected_graph_[src_num].end() &&
            undirected_graph_[dst_num].find(src_num) == undirected_graph_[dst_num].end()) {
            undirected_graph_[src_num].insert(dst_num);
            undirected_graph_[dst_num].insert(src_num);
        }
    }
    */

    /*
    void EvolutionaryTree::PrepareSubtree(std::vector<std::pair<size_t, size_t>>& edge_vector, size_t root_num) {
        if (undirected_graph_.find(root_num) != undirected_graph_.end() && !GetFlag(root_num)) {
            SetFlag(true, root_num);
            for (size_t u : undirected_graph_[root_num]) {
                //undirected_graph_[root_num].erase(u);
                //undirected_graph_[u].erase(root_num);
                if (!GetFlag(u)) {
                    edge_vector.push_back(std::make_pair(root_num, u));
                    PrepareSubtree(edge_vector, u);
                }
            }
        }
    }

    void EvolutionaryTree::PrepareSubtreeVertices(
            boost::unordered_set<size_t>& vertices_set,
            size_t root_num) {
        //we assume that component size is > 1
        if (vertices_set.find(root_num) != vertices_set.end()) {
            return;
        }
        vertices_set.insert(root_num);
        for (size_t u : undirected_graph_[root_num]) {
            PrepareSubtreeVertices(vertices_set, u);
        }
    }

    void EvolutionaryTree::PrepareSubtreeKruskal(std::vector<std::pair<size_t, size_t>>& edge_vector,
                                                 size_t root_vertex,
                                                 const annotation_utils::CDRAnnotatedCloneSet& clone_set,
                                                 EvolutionaryEdgeConstructor* edge_constructor) {
        boost::unordered_set<size_t> vertices_set;
        PrepareSubtreeVertices(vertices_set, root_vertex);
        for (size_t v : vertices_set) {
            //undirected_graph_.erase(v);
            undirected_graph_[v].clear();
        }
        size_t n = vertices_set.size();
        boost::unordered_map<size_t, size_t> vertex_to_index;
        std::vector<size_t> index_to_vertex(n);
        size_t index = 0;
        for (size_t vertex : vertices_set) {
            //std::cout << clone_to_string(clone_set[vertex]) << std::endl << std::endl;
            vertex_to_index[vertex] = index;
            index_to_vertex[index] = vertex;
            ++index;
        }
        //size_t root_index = vertex_to_index[root_vertex];
        typedef EdmondsTarjanDMSTCalculator::WeightedEdge WeightedEdge;
        std::vector<WeightedEdge> edges;
        std::vector<boost::unordered_map<size_t, size_t>> kruskal_graph(n);
        for (size_t v = 0; v < n; ++v) {
            for (size_t u = 0; u < n; ++u) {
                if (v == u) {
                    continue;
                }
                //INFO(index_to_vertex[u] << "->" << index_to_vertex[v]);
                size_t CDR3_dist = edge_constructor->ConstructEdge(clone_set[index_to_vertex[v]],
                                                                   clone_set[index_to_vertex[u]],
                                                                   index_to_vertex[v],
                                                                   index_to_vertex[u]).cdr3_distance;
                edges.push_back(WeightedEdge(v, u, static_cast<double>(CDR3_dist)));
            }
        }

        std::map<size_t, size_t> rank;
        std::map<size_t, size_t> parent;
        typedef boost::associative_property_map<std::map<size_t, size_t>> AP_map;
        boost::disjoint_sets<AP_map, AP_map> ds(
                boost::make_assoc_property_map(rank),
                boost::make_assoc_property_map(parent));
        for (size_t i = 0; i < n; ++i) {
            ds.make_set(i);
        }
        std::sort(edges.begin(), edges.end(), [](WeightedEdge e1, WeightedEdge e2) {
            return e1.weight_ < e2.weight_;
        });

        std::vector<WeightedEdge> edges_to_add;
        size_t edge_num = 0;
        size_t added_edge_num = 0;
        while (edge_num < edges.size() && added_edge_num < n) {
            auto const &edge = edges[edge_num];
            if (ds.find_set(edge.src_) == ds.find_set(edge.dst_)) {
                ++edge_num;
                continue;
            }
            edges_to_add.push_back(edge);
            ds.union_set(edge.src_, edge.dst_);
            ++edge_num;
            ++added_edge_num;
        }

        for (auto &we : edges_to_add) {
            size_t src_clone_num = index_to_vertex[we.src_];
            size_t dst_clone_num = index_to_vertex[we.dst_];
            undirected_graph_[src_clone_num].insert(dst_clone_num);
            undirected_graph_[dst_clone_num].insert(src_clone_num);
        }
        PrepareSubtree(edge_vector, root_vertex);
    }

    void EvolutionaryTree::PrepareSubtreeEdmonds(std::vector<std::pair<size_t, size_t>>& edge_vector,
                                                 size_t root_vertex,
                                                 ShmModel& model,
                                                 const annotation_utils::CDRAnnotatedCloneSet& clone_set,
                                                 EvolutionaryEdgeConstructor* edge_constructor) {
        //INFO("Preparing graph, root is " << root_vertex);
        boost::unordered_set<size_t> vertices_set;
        PrepareSubtreeVertices(vertices_set, root_vertex);
        for (size_t v : vertices_set) {
            undirected_graph_.erase(v);
        }
        size_t n = vertices_set.size();
        boost::unordered_map<size_t, size_t> vertex_to_index;

        std::vector<size_t> index_to_vertex(n);
        size_t index = 0;
        for (size_t vertex : vertices_set) {
            //std::cout << clone_to_string(clone_set[vertex]) << std::endl << std::endl;
            vertex_to_index[vertex] = index;
            index_to_vertex[index] = vertex;
            ++index;
        }
        size_t root_index = vertex_to_index[root_vertex];

        typedef EdmondsTarjanDMSTCalculator::WeightedEdge WeightedEdge;
        std::vector<WeightedEdge> edges;
        for (size_t v = 0; v < n; ++v) {
            for (size_t u = 0; u < n; ++u) {
                if (v == u) {
                    continue;
                }
                //INFO(index_to_vertex[u] << "->" << index_to_vertex[v]);
                double dist = model.CDR3TransitionProb(
                        edge_constructor->ConstructEdge(clone_set[index_to_vertex[v]],
                                                        clone_set[index_to_vertex[u]],
                                                        index_to_vertex[v],
                                                        index_to_vertex[u]));
                if (dist != 0) {
                    edges.push_back(WeightedEdge(v, u, log(dist)));
                    //std::cout << v << " " << u << " " << log(dist)  << std::endl;
                }

            }
        }
        //std::cout << "root " << root_index  << std::endl;
        //INFO("Graph prepared");
        EdmondsTarjanDMSTCalculator e_calc(n, edges, root_index);
        e_calc.EmpondsTarjan();
        std::vector<WeightedEdge> edges_to_add = e_calc.GetParentEdges();
        for (auto& we : edges_to_add) {
            size_t src_clone_num = index_to_vertex[we.src_];
            size_t dst_clone_num = index_to_vertex[we.dst_];
            edge_vector.push_back(std::make_pair(src_clone_num, dst_clone_num));
            //std::cout << we.src_ << " " << we.dst_ << " (" << we.weight_ <<  ")"  << std::endl;
        }
    }
    */

    /*
    void EvolutionaryTree::WriteEdge(const EvolutionaryEdge& edge, std::ofstream& out) { //no endl
        out << edge.src_clone->Read().id << "\t" << edge.dst_clone->Read().id << "\t"
        << edge.src_clone->Read().name << "\t" << edge.dst_clone->Read().name << "\t"
        << edge.src_clone->VSHMs().size() + edge.src_clone->JSHMs().size() << "\t"
        << edge.dst_clone->VSHMs().size() + edge.dst_clone->JSHMs().size()<< "\t"
        << edge.edge_type << "\t" <<  edge.num_intersected_shms << "\t" << edge.num_added_shms
        << "\t" << edge.cdr3_distance << "\t" << edge.weight << "\t"
        << edge.src_clone->Productive() << "\t" << edge.dst_clone->Productive() << "\t" << edge.IsSynonymous();
    }

    void EvolutionaryTree::WriteInFile(std::string output_fname) {
        std::ofstream out(output_fname);
        out << "Src_id\tDst_id\tSrc_clone\tDst_clone\tNum_Src_SHMs\tNum_Dst_SHMs\tEdge_type\t";
        out << "Num_shared_SHMs\tNum_added_SHMs\tCDR3_dist\tWeight\t";
        out << "Src_productive\tDst_productive\tSynonymous\n";
        for(auto it = edges_.begin(); it != edges_.end(); it++)
            WriteEdge(it->second, out);
        out << std::endl;
        out.close();
    }

    void EvolutionaryTree::WriteInFileWithCDR3s(std::string output_fname) {
        std::ofstream out(output_fname);
        out << "Src_id\tDst_id\tSrc_clone\tDst_clone\tNum_Src_SHMs\tNum_Dst_SHMs\tEdge_type\t";
        out << "Num_shared_SHMs\tNum_added_SHMs\tCDR3_dist\tWeight\t";
        out << "Src_productive\tDst_productive\tSynonymous\t";
        out << "Src_CDR3\tDst_CDR3\n";
        for(auto it = edges_.begin(); it != edges_.end(); it++) {
            WriteEdge(it->second, out);
            std::string src_CDR3_string;
            std::string dst_CDR3_string;
            auto const& src_CDR3 = it->second.src_clone->CDR3();
            auto const& dst_CDR3 = it->second.dst_clone->CDR3();
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

    void EvolutionaryTree::WriteVerticesInFile(std::string output_fname,
                                               const annotation_utils::CDRAnnotatedCloneSet& clone_set) {
        std::ofstream out(output_fname);
        out << "Clone_id\tClone_name\tProductive\tAA_seq\tOFR\tLeft_CDR3_anchor_AA\tRight_CDR3_anchor_AA\n";
        for (auto p : passed_flag_) {
            auto const& clone = clone_set[p.first];
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
    */

    bool EvolutionaryTree::IsRoot(size_t clone_id) const {
        VERIFY_MSG(vertices_.find(clone_id) != vertices_.end(), "Tree does not contain vertex " << clone_id);
        return edges_.find(clone_id) == edges_.end();
    }

    size_t EvolutionaryTree::GetRoot() const {
        for(auto it = vertices_.begin(); it != vertices_.end(); it++) {
            if(IsRoot(*it)) {
                return *it;
            }
        }
        VERIFY_MSG(false, "Root was not found");
        return size_t(-1);
    }

    bool EvolutionaryTree::IsLeaf(size_t clone_id) const {
        VERIFY_MSG(ContainsClone(clone_id), "Tree does not contain vertex " << clone_id);
        return outgoing_edges_.find(clone_id) == outgoing_edges_.end();
    }

    bool EvolutionaryTree::IsForest() const {
        size_t num_roots = 0;
        for(auto it = vertices_.begin(); it != vertices_.end(); it++) {
            if(IsRoot(*it)) {
                num_roots++;
            }
        }
        return num_roots > 1;
    }

    size_t EvolutionaryTree::GetRootNumber() const {
        size_t num_roots = 0;
        for(auto it = vertices_.begin(); it != vertices_.end(); it++) {
            if(IsRoot(*it)) {
                num_roots++;
            }
        }
        return num_roots;
    }

    size_t EvolutionaryTree::EdgeDepth() const {
        VERIFY_MSG(false, "Implement me!");
        return 0;
    }

    const std::vector<EvolutionaryEdge>& EvolutionaryTree::OutgoingEdges(size_t clone_id) const {
        VERIFY_MSG(ContainsClone(clone_id), "Tree does not contain vertex " << clone_id);
        VERIFY_MSG(!IsLeaf(clone_id), "Vertex " << clone_id << " is leaf");
        return outgoing_edges_.at(clone_id);
    }

    std::vector<size_t> EvolutionaryTree::GetRoots() const {
        std::vector<size_t> roots;
        for(auto it = vertices_.begin(); it != vertices_.end(); it++) {
            if(IsRoot(*it)) {
                roots.push_back(*it);
            }
        }
        return roots;
    }

    void EvolutionaryTree::SetTreeOutputFname(std::string output_dir, size_t index1, size_t index2, size_t v_num, size_t e_num) {
        std::stringstream ss;
        ss << "clonal_tree_" << index1 << "-" << index2 << "_Vsize_" << v_num << "_Esize_" << e_num << ".tree";
        tree_output_fname_ = path::append_path(output_dir, ss.str());
    }
    std::string EvolutionaryTree::GetTreeOutputFname() const{
        return tree_output_fname_;
    }
    void EvolutionaryTree::SetVerticesOutputFname(std::string output_dir, size_t index1, size_t index2, size_t v_num, size_t e_num) {
        std::stringstream ss;
        ss << "clonal_tree_" << index1 << "-" << index2 << "_Vsize_" << v_num << "_Esize_" << e_num << ".tree";
        vertices_output_fname_ = path::append_path(output_dir, ss.str());
    }
    std::string EvolutionaryTree::GetVerticesOutputFname() const{
        return vertices_output_fname_;
    }


    std::ostream& operator<<(std::ostream& out, const EvolutionaryTree &tree) {
        for(auto it = tree.cbegin(); it != tree.cend(); it++)
            out << *it << std::endl;
        return out;
    }
}
