#pragma once

#include <memory>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include "spliced_read.hpp"
#include <iostream>

#define DISTANCE_THRESHOLD 2

namespace ig_repertoire_constructor{

struct Node;
typedef std::shared_ptr <Node> NodePtr;

struct Node {
    Node(size_t read_ind) : child1(nullptr), child2(nullptr), read_ind(read_ind), size(1) { }
    Node(NodePtr child1, NodePtr child2, double dist) :
        child1(child1), child2(child2), read_ind(0), size(child1->size + child2->size), distance(dist) { }

    void UpdateElements(std::vector <size_t> &elements) {
        if (child1) {
            child1->UpdateElements(elements);
            child2->UpdateElements(elements);
        } else {
            elements.push_back(read_ind);
        }
    }

    std::string str() {
        if (child1) {
            return "*" + std::to_string(distance) + "(" + child1->str() + " " + child2->str() + ")";
        } else {
            return std::to_string(read_ind);
        }
    }

    NodePtr child1;
    NodePtr child2;
    size_t read_ind;
    int size;
    double distance;
};

struct Hash{
    size_t operator()(const std::pair <NodePtr, NodePtr> &key) const {
        return std::hash<NodePtr>()(key.first) + std::hash<NodePtr>()(key.second);
    }
};

class UPGMAClusterization {
public:
    UPGMAClusterization() { }

    std::shared_ptr <std::vector< std::vector <size_t> > > Clusterize(std::vector <SplicedRead> const &spliced_reads);

private:
    struct QueueElement {
        bool operator <(QueueElement const & e) const {
            return distance > e.distance;
        }

        NodePtr node1;
        NodePtr node2;
        double distance;
    };

    typedef std::pair <NodePtr, double> node_distance_pair;
    void InitNodes(std::vector <SplicedRead> const &spliced_reads);
    void CountDistances(std::vector <SplicedRead> const &spliced_reads);
    std::pair <NodePtr, NodePtr> GetClosestClusters();
    void MergeAndUpdateDistances(NodePtr i, NodePtr j);

    std::unordered_set <NodePtr> nodes_;
    std::priority_queue <QueueElement> distances_heap_;
    std::unordered_map <std::pair <NodePtr, NodePtr>, double, Hash> distances_;
};

std::shared_ptr <std::vector <std::vector <size_t> > > UPGMAClusterization::Clusterize(std::vector <SplicedRead> const &spliced_reads) {
    InitNodes(spliced_reads);
    CountDistances(spliced_reads);
    while (nodes_.size() > 1) {
        std::pair <NodePtr, NodePtr> closest = GetClosestClusters();
        MergeAndUpdateDistances(closest.first, closest.second);
    }

    /*
    if (spliced_reads.size() > 3) {
        INFO(nodes_[0]->str());
    }
    */

    std::vector <NodePtr> cluster_nodes;
    std::vector <NodePtr> queue;
    queue.push_back(*(nodes_.begin()));
    while (!queue.empty()) {
        NodePtr curr = queue.back();
        queue.pop_back();
        if (curr->distance < DISTANCE_THRESHOLD || !curr->child1 || !curr->child1->child1 || !curr->child2->child1) {
            cluster_nodes.push_back(curr);
        } else {
            queue.push_back(curr->child1);
            queue.push_back(curr->child2);
        }
    }

    std::shared_ptr <std::vector <std::vector <size_t> > > clusters_ptr(new std::vector <std::vector <size_t> >(cluster_nodes.size()));
    for (size_t i = 0; i != cluster_nodes.size(); ++i) {
        cluster_nodes[i]->UpdateElements((*clusters_ptr)[i]);
    }
    return clusters_ptr;
}

void UPGMAClusterization::InitNodes(std::vector<SplicedRead> const &spliced_reads) {
    for (size_t i = 0; i != spliced_reads.size(); ++i) {
        nodes_.insert(NodePtr(new Node(i)));
    }
}

void UPGMAClusterization::CountDistances(std::vector<SplicedRead> const &spliced_reads) {
    int min_pos = GetLeftmostPosition(spliced_reads, std::vector <size_t>());
    int max_pos = GetRightmostPosition(spliced_reads, std::vector <size_t>());
    std::vector <std::vector <int> > counts(max_pos - min_pos + 1, {0, 0, 0, 0});
    bool cont = true;
    for (int pos = min_pos; cont; ++pos) {
        cont = false;
        for (size_t i = 0; i != spliced_reads.size(); ++i) {
            if (spliced_reads[i].HasCharWithIndex(pos)) {
                cont = true;
                ++counts[pos - min_pos][spliced_reads[i][pos]];
            }
        }
    }

    for (auto it = nodes_.begin(); it != nodes_.end(); ++it) {
        for (auto jt = it; jt != nodes_.end(); ++jt) {
            if (it == jt) {
                continue;
            }
            double dist = 0;
            cont = true;
            for (int pos = min_pos; cont; ++pos) {
                if (!spliced_reads[(*it)->read_ind].HasCharWithIndex(pos) || !spliced_reads[(*jt)->read_ind].HasCharWithIndex(pos)) {
                    if (pos > 0) {
                        cont = false;
                    }
                    continue;
                }
                char ci = spliced_reads[(*it)->read_ind][pos];
                char cj = spliced_reads[(*jt)->read_ind][pos];
                if (ci != cj && counts[pos - min_pos][ci] != 1 && counts[pos - min_pos][cj] != 1) {
                    ++dist;
                }
            }
            if (*it < *jt) {
                distances_[std::make_pair(*it, *jt)] = dist;
            } else {
                distances_[std::make_pair(*jt, *it)] = dist;
            }
            distances_heap_.push({*it, *jt, dist});
        }
    }
}

std::pair<NodePtr, NodePtr> UPGMAClusterization::GetClosestClusters() {
    NodePtr n1(nullptr);
    NodePtr n2(nullptr);
    while (!distances_heap_.empty() && (nodes_.find(n1) == nodes_.end() || nodes_.find(n2) == nodes_.end())) {
        n1 = distances_heap_.top().node1;
        n2 = distances_heap_.top().node2;
        distances_heap_.pop();
    }
    return std::make_pair(n1, n2);
}

void UPGMAClusterization::MergeAndUpdateDistances(NodePtr i, NodePtr j) {
    auto dist_ij_it = (i < j) ? distances_.find(std::make_pair(i, j)) :
                                distances_.find(std::make_pair(j, i));
    NodePtr new_node(new Node(i, j, dist_ij_it->second));
    distances_.erase(dist_ij_it);
    nodes_.erase(i);
    nodes_.erase(j);
    for (auto node_it = nodes_.begin(); node_it != nodes_.end(); ++node_it) {
        auto dist_i_it = (i < *node_it) ? distances_.find(std::make_pair(i, *node_it)) :
                                          distances_.find(std::make_pair(*node_it, i));
        auto dist_j_it = (j < *node_it) ? distances_.find(std::make_pair(j, *node_it)) :
                                          distances_.find(std::make_pair(*node_it, j));
        double dist = (i->size * dist_i_it->second + j->size * dist_j_it->second) / (i->size + j->size);

        distances_heap_.push({*node_it, new_node, dist});
        if (new_node < *node_it) {
            distances_[std::make_pair(new_node, *node_it)] = dist;
        } else {
            distances_[std::make_pair(*node_it, new_node)] = dist;
        }
        distances_.erase(dist_i_it);
        distances_.erase(dist_j_it);
    }
    nodes_.insert(new_node);
}

}
