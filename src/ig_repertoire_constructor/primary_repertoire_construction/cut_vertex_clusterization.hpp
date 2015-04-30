#pragma once

#include <unordered_set>
#include <unordered_map>
#include <sstream>
#include "spliced_read.hpp"

namespace ig_repertoire_constructor {

struct Edge {
    Edge(size_t from_ind, size_t to_ind, size_t weight, size_t num_patterns) :
            from_ind(from_ind),
            to_ind(to_ind),
            weight(weight),
            num_patterns(num_patterns),
            color(0) { }

    size_t from_ind;
    size_t to_ind;
    size_t weight;
    size_t num_patterns;

    size_t color;
};

struct EdgeWeight {
    size_t dist;
    size_t num_patterns;
    vector<size_t> shm_positions;

    EdgeWeight(size_t new_dist, size_t new_num_patterns, vector<size_t> shm_positions):
            dist(new_dist),
            num_patterns(new_num_patterns),
            shm_positions(shm_positions) { }
};

typedef std::shared_ptr <Edge> EdgePtr;

class CutVertexClusterization {
public:
    CutVertexClusterization(size_t max_dist_) :
            max_distance_(max_dist_),
            clusters_number_(0),
            timer_(0) { }

    std::shared_ptr <std::vector <std::vector <size_t> > > Clusterize(
            std::vector <SplicedRead> const &spliced_reads, size_t group_id);

private:
    struct AdjClusterInfo {
        AdjClusterInfo(size_t degree, size_t min_distance)
            :degree(degree), min_distance(min_distance) { }
        AdjClusterInfo()
            :degree(0), min_distance(1000000) { }
        bool operator<(AdjClusterInfo const &a) const {
            if (degree != a.degree) {
                return degree < a.degree;
            }
            return min_distance > a.min_distance;
        }

        size_t degree;
        size_t min_distance;
    };

    struct IsEmptyVectorPredicate {
        bool operator()(std::vector <size_t> const &v) const {
            return v.empty();
        }
    };

    struct DFSStackElement {
        size_t v, p;
        EdgePtr edge_from_parent;
        std::vector <EdgePtr>::iterator edge_it;
    };

    void BuildHammingGraph(std::vector <SplicedRead> const &spliced_reads);
    void PrintHammingGraph(size_t group_id) const;
    EdgeWeight GetHammingDistance(SplicedRead const &r1, SplicedRead const &r2) const;
    void DFSNonRecursive(size_t v);
    void StartDFS(size_t v, size_t p, EdgePtr edge_from_parent, std::stack<DFSStackElement> &dfs_stack);
    void EndDFS(size_t v, size_t p, EdgePtr edge_from_parent);
    void ProcessEdge(size_t v, size_t p, EdgePtr curr_edge, std::stack<DFSStackElement> &dfs_stack);

    // void DFS(size_t v, size_t p);
    void PaintEdges();
    void AssignColorToVectex();

    size_t max_distance_;
    size_t clusters_number_;
    std::vector <std::vector <EdgePtr> > graph_;
    size_t timer_;
    std::vector <size_t> used_;
    std::vector <size_t> tenter_;
    std::vector <size_t> treturn_;
    std::stack <EdgePtr, std::vector <EdgePtr> > edges_stack_;
    std::vector <size_t> vertex_colors_;
};

std::shared_ptr<std::vector<std::vector<size_t> > > CutVertexClusterization::Clusterize(
            const std::vector<SplicedRead> &spliced_reads, size_t group_id) {
    if (spliced_reads.size() < 4) {
        std::shared_ptr <std::vector <std::vector <size_t > > > res(new std::vector <std::vector <size_t> >(1, std::vector <size_t>()));
        for (size_t i = 0; i != spliced_reads.size(); ++i) {
            (*res)[0].push_back(i);
        }
        return res;
    }
    BuildHammingGraph(spliced_reads);
    PrintHammingGraph(group_id);
    PaintEdges();

    /*
    std::vector <size_t> colors;
    for (size_t v = 0; v != graph_.size(); ++v) {
        for (size_t e = 0; e != graph_[v].size(); ++e) {
            colors.push_back(graph_[v][e]->color);
        }
    }
    INFO(colors);
    */

    AssignColorToVectex();
    std::shared_ptr <std::vector <std::vector <size_t > > > res(new std::vector <std::vector <size_t> >(clusters_number_, std::vector <size_t>()));
    for (size_t v = 0; v != graph_.size(); ++v) {
        size_t curr_color = vertex_colors_[v];
        if (graph_[v].size() == 1) {
            size_t to = (graph_[v][0]->from_ind == v) ? graph_[v][0]->to_ind : graph_[v][0]->from_ind;
            curr_color = vertex_colors_[to];
        }
        (*res)[curr_color].push_back(v);
    }

    res->erase(std::remove_if(res->begin(), res->end(), IsEmptyVectorPredicate()), res->end());

    return res;
}

void CutVertexClusterization::BuildHammingGraph(const std::vector<SplicedRead> &spliced_reads) {
    graph_.assign(spliced_reads.size(), std::vector <EdgePtr>());

    // tmp
    vector<set<size_t> > shm_positions;
    vector<size_t> reads_ids;
    for(size_t i = 0; i < spliced_reads.size(); i++) {
        reads_ids.push_back(i);
        shm_positions.push_back(set<size_t>());
    }

    for (size_t i = 0; i != spliced_reads.size(); ++i) {
        for (size_t j = i + 1; j != spliced_reads.size(); ++j) {
            auto weight = GetHammingDistance(spliced_reads[i], spliced_reads[j]);
            if (weight.dist <= max_distance_) {
                EdgePtr new_edge(new Edge(i, j, weight.dist, weight.num_patterns));
                graph_[i].push_back(new_edge);
                graph_[j].push_back(new_edge);

                // tmp
                for(auto it = weight.shm_positions.begin(); it != weight.shm_positions.end(); it++) {
                    shm_positions[i].insert(*it);
                    shm_positions[j].insert(*it);
                }

                if(weight.dist == 0)
                    reads_ids[j] = i;
            }
        }
    }

    // write in file
    stringstream ss;
    ss << "shm_positions_tau" << ig_cfg::get().aligning_params.overlap_mismatches_threshold << ".txt";
    string fname = ss.str();
    ofstream out(fname, std::ios_base::app);
    for(size_t i = 0; i < shm_positions.size(); i++) {
        if(reads_ids[i] != i)
            continue;
        out << spliced_reads[i].ReadName() << "\t" << spliced_reads[i].ReadLength() << endl;
        for(auto it = shm_positions[i].begin(); it != shm_positions[i].end(); it++)
            out << spliced_reads[i].AbsoletePosition(*it) << "\t";
        out << endl;
    }

    out.close();
}

string GetHgraphName(size_t group_id, size_t graph_size) {
    stringstream ss;
    ss << "hgraph_" << group_id << "_size_" << graph_size;
    return ss.str();
}

void CutVertexClusterization::PrintHammingGraph(size_t group_id) const {
    if (!ig_cfg::get().io.output_ham_graphs) {
        return;
    }

    size_t min_size = 100;
    if(graph_.size() < min_size)
        return;

    string output_fname = GetHgraphName(group_id, graph_.size());
    std::string file_name = path::append_path(ig_cfg::get().io.hgraph_dir, output_fname);
    INFO("Saving Hamming graph into '" << file_name << "'");
    std::fstream fs(file_name.c_str(), std::ios::out);
    VERIFY(fs.good());
    fs << "n=" << graph_.size() << "\n";
    for (size_t i = 0; i != graph_.size(); ++i) {
        for (const EdgePtr & edge : graph_[i]) {
            if (i == edge->from_ind) {
                fs << edge->from_ind << " " <<  edge->to_ind << " " << edge->weight << " " <<
                        edge->num_patterns << "\n";
            }
        }
    }
    VERIFY(!fs.fail());
    fs.close();
}


bool PositionIsPattern(const SplicedRead &r, size_t i, char N1, char N2) {
    return nucl(r[i]) == N1 or nucl(r[i]) == N2;
}

bool PositionIsRGYM(const SplicedRead &r, size_t i) {
    //cout << "-------------------" << endl;
    //cout << "RGYM check starts" << endl;
    // i is border position
    if(i == 0 or i >= r.size() - 2) {
    //    cout << "RGYM check. Wrong index" << endl;
        return false;
    }
    //cout << "RGYM: " << nucl(r[i - 1]) << nucl(r[i]) << nucl(r[i + 1]) << nucl(r[i + 2]) << endl;

    // hot spot is G
    if(nucl(r[i]) != 'G') {
    //    cout << "RGYM. Hotspot is not G" << endl;
        return false;
    }

    // A-T at i - 1
    // C-T at i + 1
    // A-T at i + 2
    bool res = PositionIsPattern(r, i - 1, 'A', 'G') and
            PositionIsPattern(r, i + 1, 'C', 'T') and
            PositionIsPattern(r, i + 2, 'A', 'T');

    //cout << "RGYM: A/G G C/T A/T: " << res << endl;
    return res;
}

bool PositionIsWRCY(const SplicedRead &r, size_t i) {
    //cout << "-------------------" << endl;
    //cout << "WRCY check starts" << endl;
    if(i <= 2 or i == r.size() - 1) {
    //    cout << "WRCY. Wrong index" << endl;
        return false;
    }
    //cout << "WRCY: " << nucl(r[i - 2]) << nucl(r[i - 1]) << nucl(r[i]) << nucl(r[i + 1]) << endl;

    // hot spot is C
    if(nucl(r[i]) != 'C') {
    //    cout << "WRCY. Hotspot is not C" << endl;
        return false;
    }

    // A-T at i - 2
    // A-G at i - 1
    // C-T at i + 1
    bool res = PositionIsPattern(r, i - 2, 'A', 'T') and
            PositionIsPattern(r, i - 1, 'A', 'G') and
            PositionIsPattern(r, i + 1, 'C', 'T');

    //cout << "WRCY: A/T A/G C C/T: " << res << endl;
    return res;
}

bool IsMismatchSHMPattern(const SplicedRead &r1, const SplicedRead &r2, size_t i) {
    // position is not mismatch
    if(r1[i] == r2[i])
        return false;

    if(PositionIsWRCY(r1, i) or PositionIsWRCY(r2, i))
        return true;

    return PositionIsRGYM(r1, i) or PositionIsRGYM(r2, i);
}

EdgeWeight CutVertexClusterization::GetHammingDistance(const SplicedRead &r1, const SplicedRead &r2) const {
    size_t dist = 0;
    size_t num_patterns = 0;
    vector<size_t> shm_positions;
    for (size_t i = 0; i != r1.size(); ++i) {
        if (r1[i] != r2[i]) {
            ++dist;
            if(IsMismatchSHMPattern(r1, r2, i)) {
                num_patterns++;
                shm_positions.push_back(i);
            }
            if (dist > max_distance_) {
                return EdgeWeight(dist, num_patterns, shm_positions);
            }
        }
    }
    return EdgeWeight(dist, num_patterns, shm_positions);
}


void CutVertexClusterization::DFSNonRecursive(size_t start) {
    std::stack <DFSStackElement> my_dfs_stack;
    StartDFS(start, graph_.size(), nullptr, my_dfs_stack);
    while (!my_dfs_stack.empty()) {
        DFSStackElement & top = my_dfs_stack.top();
        size_t v = top.v;
        size_t p = top.p;
        if (top.edge_it == graph_[v].end()) {
            EndDFS(v, p, top.edge_from_parent);
            my_dfs_stack.pop();
        } else {
            EdgePtr curr_edge = *top.edge_it;
            top.edge_it++;
            ProcessEdge(v, p, curr_edge, my_dfs_stack);
        }
    }
}

void CutVertexClusterization::StartDFS(size_t v, size_t p, EdgePtr edge_from_parent, std::stack <DFSStackElement> &dfs_stack) {
    used_[v] = true;
    ++timer_;
    tenter_[v] = timer_;
    treturn_[v] = timer_;
    dfs_stack.push({v, p, edge_from_parent, graph_[v].begin()});
}

void CutVertexClusterization::EndDFS(size_t v, size_t p, EdgePtr edge_from_parent) {
    if (!edge_from_parent) {
        return;
    }
    if (treturn_[v] >= tenter_[p]) {
        size_t color = clusters_number_;
        ++clusters_number_;
        while (edges_stack_.top() != edge_from_parent) {
            edges_stack_.top()->color = color;
            edges_stack_.pop();
        }
        edge_from_parent->color = color;
        edges_stack_.pop();
    }
    if (treturn_[v] < treturn_[p]) {
        treturn_[p] = treturn_[v];
    }

}

void CutVertexClusterization::ProcessEdge(size_t v, size_t p, EdgePtr curr_edge, std::stack <DFSStackElement> &dfs_stack){
    size_t to = (curr_edge->from_ind == v) ? curr_edge->to_ind : curr_edge->from_ind;
    if (to == p) {
        return;
    }
    if (!used_[to]) {
        edges_stack_.push(curr_edge);
        StartDFS(to, v, curr_edge, dfs_stack);
   } else {
        if (tenter_[to] < tenter_[v]) {
            edges_stack_.push(curr_edge);
        }
        if (treturn_[v] > tenter_[to]){
            treturn_[v] = treturn_[to];
        }
    }

}

/*
void CutVertexClusterization::DFS(size_t v, size_t p) {
    used_[v] = true;
    ++timer_;
    tenter_[v] = timer_;
    treturn_[v] = timer_;
    for (auto edge_it = graph_[v].begin(); edge_it != graph_[v].end(); ++edge_it) {
        EdgePtr curr_edge = *edge_it;
        size_t to = (curr_edge->from_ind == v) ? curr_edge->to_ind : curr_edge->from_ind;
        if (to == p) {
            continue;
        }
        if (!used_[to]) {
            edges_stack_.push(curr_edge);
            DFS(to, v);
            if (treturn_[to] >= tenter_[v]) {
                size_t color = clusters_number_;
                ++clusters_number_;
                while (edges_stack_.top() != curr_edge) {
                    edges_stack_.top()->color = color;
                    edges_stack_.pop();
                }
                curr_edge->color = color;
                edges_stack_.pop();
            }
            if (treturn_[to] < treturn_[v]) {
                treturn_[v] = treturn_[to];
            }
        } else {
            if (tenter_[to] < tenter_[v]) {
                edges_stack_.push(curr_edge);
            }
            if (treturn_[v] > tenter_[to]){
                treturn_[v] = treturn_[to];
            }
        }
    }
}
*/

void CutVertexClusterization::PaintEdges() {
    used_.assign(graph_.size(), false);
    treturn_.assign(graph_.size(), 0);
    tenter_.assign(graph_.size(), 0);
    for (size_t v = 0; v != graph_.size(); ++v) {
        if (!used_[v]) {
            timer_ = 0;
            // DFS(v, graph_.size());
            DFSNonRecursive(v);
        }
    }
}

void CutVertexClusterization::AssignColorToVectex() {
    vertex_colors_.assign(graph_.size(), 0);
    for (size_t v = 0; v != graph_.size(); ++v) {
        std::unordered_map <size_t, AdjClusterInfo> cluster_adj;
        for (auto edge_it = graph_[v].begin(); edge_it != graph_[v].end(); ++edge_it) {
            EdgePtr curr_edge = *edge_it;
            size_t curr_color_id = curr_edge->color;
            if (cluster_adj.find(curr_color_id) == cluster_adj.end()) {
                cluster_adj[curr_color_id] = AdjClusterInfo(0, max_distance_ + 1);
            }
            auto curr_adj_it = cluster_adj.find(curr_color_id);
            ++(curr_adj_it->second.degree);
            curr_adj_it->second.min_distance = std::min(curr_adj_it->second.min_distance, (size_t)curr_edge->weight);
        }

        auto best_cluster_it = cluster_adj.begin();
        for (auto adj_cluster_it = cluster_adj.begin(); adj_cluster_it != cluster_adj.end(); ++adj_cluster_it) {
            if (best_cluster_it->second < adj_cluster_it->second) {
                best_cluster_it = adj_cluster_it;
            }
        }
        vertex_colors_[v] = best_cluster_it->first;
    }
}

}
