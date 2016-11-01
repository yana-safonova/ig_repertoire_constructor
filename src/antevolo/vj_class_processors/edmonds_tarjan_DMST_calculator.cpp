#include <logger/logger.hpp>
#include "edmonds_tarjan_DMST_calculator.hpp"

namespace antevolo {
    void EdmondsTarjanDMSTCalculator::InitVertex() {
        ++vertices_num_; // the nubmer of vertices
        in_.push_back(WeightedEdge());
        const_add_.push_back(-1);
        parent_.push_back(size_t(-1));
        phase_.push_back(size_t(-1));
        children_.push_back(std::vector<size_t>());
        Ps_.push_back(P_queue());
    }

    void EdmondsTarjanDMSTCalculator::Initialize() {
        for (size_t u = 0; u < n; ++u) {
            InitVertex();
        }
        for (auto we : edges_) {
            if (we.dst_ != root_) {
                Ps_[we.dst_].push(we);
                edge_initial_weight_map_[std::make_pair(we.src_, we.dst_)] = we.weight_;
            }
        }
        TRACE("Initialization finished");
    }
    /*
    void EdmondsTarjanDMSTCalculator::Contract(boost::disjoint_sets<AP_map, AP_map> &ds_contract,
                                               boost::disjoint_sets<AP_map, AP_map> &ds_arbor) {
        INFO("Contraction phase starts");
        Initialize();
        for (auto p : Ps_) {
            if (p.empty()) continue;
            std::cout << p.top().weight_ << " ";
        }  std::cout << std::endl;
        for (size_t a_start = 0; a_start < n; ++a_start) {
            size_t a = a_start;
            INFO("vertex " << a);
            while (!Ps_[a].empty()) {
                INFO(a << " p_queue size = " << Ps_[a].size());
                WeightedEdge e = Ps_[a].top();
                //INFO("edge " << e.src_ << "->" << e.dst_ << " (" << e.weight_ << ")");
                size_t b = parent_label_map_[ds_contract.find_set(e.src_)];
                //INFO("src parent vertex " << b);
                if (a == b) { //if e is a selfloop
                    Ps_[a].pop();
                    continue;
                }
                if (in_[a] == e) { //if e is an already set edge
                    break;
                }
                in_[a] = e;
                prev_[a] = b;
                //INFO("src in_edge " << in_[e.src_].src_ << "->" << in_[e.src_].dst_ << " (" << in_[e.src_].weight_ << ")");
                if (ds_arbor.find_set(a) != ds_arbor.find_set(e.src_)) {
                    ds_arbor.union_set(a, e.src_);
                    a = b;
                    continue;
                }
                InitVertex();
                size_t c = vertices_num_-1;
                ds_contract.make_set(c);
                ds_arbor.make_set(c);
                std::vector<size_t> cycle;
                cycle.push_back(a);
                INFO(a);
                while (parent_label_map_[ds_contract.find_set(in_[a].src_)] != cycle[0]) {
                    INFO(c << " | " << a << " " << parent_label_map_[ds_contract.find_set(in_[a].src_)]);
                    a = parent_label_map_[ds_contract.find_set(in_[a].src_)]; // HERE
                    cycle.push_back(a);
                }
                for (auto a2 : cycle) {
                    INFO(c << " || " << a2 << " size " << Ps_[a2].size());
                    parent_[a2] = c;
                    const_add_[a2] = -in_[e.src_].weight_; // - ?
                    //INFO(a2 << " const add " << const_add_[a2]);
                    children_[c].push_back(a2);
                    ds_contract.union_set(c, a2);
                    ds_arbor.union_set(c, a2);
                    parent_label_map_[ds_contract.find_set(c)] = c;
                    while (!Ps_[a2].empty()) {
                        INFO(c << " " << a2 << " const add " << in_[a2].weight_);
                        WeightedEdge e2 = Ps_[a2].top();
                        //if (a2 == e.) {
                        //    e2.weight_ -= in_[a2].weight_;
                        //}
                        e2.weight_ -= in_[a2].weight_;
                        Ps_[a2].pop();
                        Ps_[c].push(e2);
                    }
                    //a = prev_[a];
                }
                //e.weight_ -= in_[a]
                //Ps_[c].push(e);
                WeightedEdge we = Ps_[c].top();
                INFO("edge " << we.src_ << "->" << we.dst_ << " (" << we.weight_ << ")");
                a = c;// sure?
            }
        }
        INFO("Contraction finished")
        for (size_t v = 0; v < children_.size(); ++v) {
            std::cout << v << " | ";
            for (auto u : children_[v]) {
                std::cout << u << " ";
            }
            std::cout << std::endl;
        }
    }

     */

    void EdmondsTarjanDMSTCalculator::Contract2(boost::disjoint_sets<AP_map, AP_map> &ds_contract) {
        TRACE("Contraction phase starts");
        Initialize();
        /*
        for (auto p : Ps_) {
            if (p.empty()) continue;
            std::cout << p.top().weight_ << " ";
        }  std::cout << std::endl;
        */
        //std::vector<size_t> next(vertices_num_, size_t(-1));
        phase_[root_] = 0;
        for (size_t a_start = 0; a_start < n; ++a_start) {
            size_t crnt_phase = a_start+1;
            size_t a = parent_label_map_[ds_contract.find_set(a_start)];
            if (phase_[a] != size_t(-1)) {
                continue;
            }
            //INFO("vertex " << a);
            //next.assign(vertices_num_, size_t(-1));
            while (!Ps_[a].empty()) {
                /*
                for (size_t i = 0; i < vertices_num_; ++i) {
                    std::cout << parent_label_map_[ds_contract.find_set(i)] << " ";
                } std::cout << std::endl;
                for (size_t i = 0; i < vertices_num_; ++i) {
                    if (next[i] != size_t(-1)) {
                        std::cout << i << " " << next[i] << " | ";
                    }
                } std::cout << std::endl;
               */
                phase_[a] = crnt_phase;

                //INFO(a << " p_queue size = " << Ps_[a].size());
                WeightedEdge e = Ps_[a].top();
                size_t b = parent_label_map_[ds_contract.find_set(e.src_)];
                if (a == b) { //if e is a selfloop
                    Ps_[a].pop();
                    continue;
                }
                //next[a] = b;
                //INFO("next[" << a << "] = " << b);
                in_[a] = e;
                //INFO("src in_edge " << in_[e.src_].src_ << "->" << in_[e.src_].dst_ << " (" << in_[e.src_].weight_ << ")");

                //if (next[b] == size_t(-1)) {
                if (phase_[b] == size_t(-1)) {
                    a = b;
                    continue;
                }
                if (phase_[b] != crnt_phase) { //if b has been already rooted
                    break;
                }

                InitVertex();
                //next.push_back(size_t(-1));
                size_t c = vertices_num_-1;
                ds_contract.make_set(c);
                std::vector<size_t> cycle;
                cycle.push_back(a);
                //INFO(a);
                size_t u = parent_label_map_[ds_contract.find_set(in_[a].src_)];
                while (u != cycle[0]) {
                    //INFO(c << " | " << u << " " << parent_label_map_[ds_contract.find_set(next[u])]);
                    cycle.push_back(u);
                    u = parent_label_map_[ds_contract.find_set(in_[u].src_)];
                }
                for (auto a2 : cycle) {
                    ds_contract.union_set(c, a2);
                    parent_label_map_[ds_contract.find_set(c)] = c;
                    //INFO(c << " || " << a2 << " size " << Ps_[a2].size());
                    parent_[a2] = c;
                    //INFO(a2 << " const add " << const_add_[a2]);
                    children_[c].push_back(a2);
                    while (!Ps_[a2].empty()) {
                        //INFO(c << " " << a2 << " const add " << in_[a2].weight_);
                        WeightedEdge e2 = Ps_[a2].top();
                        if (parent_label_map_[ds_contract.find_set(e2.dst_)] != c) {
                            INFO("invalid in_edge, top level: " << c << " bottom level " << a2 << " src " << e2.src_
                                 << " dst " << e2.dst_);
                        }
                        if (parent_label_map_[ds_contract.find_set(e2.src_)] == c) {
                            Ps_[a2].pop();
                            continue;
                        }
                        //if (a2 == e.) {
                        //    e2.weight_ -= in_[a2].weight_;
                        //}
                        e2.weight_ -= in_[a2].weight_;
                        Ps_[a2].pop();
                        Ps_[c].push(e2);


                    }
                    //a = prev_[a];
                }

                //WeightedEdge we = Ps_[c].top();
                //INFO("edge " << we.src_ << "->" << we.dst_ << " (" << we.weight_ << ")");
                a = c;// sure?
            }
        }
        TRACE("Contraction finished")
        /*
        for (size_t v = 0; v < children_.size(); ++v) {
            std::cout << v << " | ";
            for (auto u : children_[v]) {
                std::cout << u << " ";
            }
            std::cout << " in_edge " << in_[v].src_ << "->" << in_[v].dst_ << " (" << in_[v].weight_ << ")" << std::endl;
        }
        */
    }
    /*
    void EdmondsTarjanDMSTCalculator::Dismantle(size_t u) {
        //INFO("dismantling " << u);
        while (parent_[u] != size_t(-1)) {
            u = parent_[u];
            for (size_t v : children_[u]) {
                parent_[v] = size_t(-1);
                if (children_[v].size() != 0) {
                    R_.push(v);
                }
            }
            children_[u].clear();
        }
    }

    void EdmondsTarjanDMSTCalculator::Expand(size_t r) {
        R_ = std::stack<size_t>();
        Dismantle(r);
        while (R_.size() != 0) {
            size_t c = R_.top();
            R_.pop();
            WeightedEdge e = in_[c];
            in_[e.dst_] = e;
            Dismantle(e.dst_);
        }
    }
    */

    void EdmondsTarjanDMSTCalculator::InEdgeHandle(size_t u) {
        //INFO("InEdgeHandling " << u);
        //WeightedEdge in_edge = in_[u];
        //INFO("edge " << in_edge.src_ << "->" << in_edge.dst_ << " (" << in_edge.weight_ << ")");
        //std::vector<size_t> path;
        size_t v = in_[u].dst_;
        while (parent_[v] != u) {
            v = parent_[v];
        }
        in_[v] = in_[u];
    }


    void EdmondsTarjanDMSTCalculator::ExpandWithEdgeHandling()
    {
        for (size_t r = vertices_num_-1; r >= n; --r) {
            InEdgeHandle(r);
        }
    }

    void EdmondsTarjanDMSTCalculator::EmpondsTarjan() {
        std::map<size_t, size_t> rank_contract;
        std::map<size_t, size_t> parent_contract;
        // boost::disjoint_sets<AP_map, AP_map> dson_undirected_edges(rank, parent);
        boost::disjoint_sets<AP_map, AP_map> ds_contract(
                boost::make_assoc_property_map(rank_contract),
                boost::make_assoc_property_map(parent_contract));
        for (size_t i = 0; i < n; ++i) {
            ds_contract.make_set(i);
            parent_label_map_[i] = i;
        }

        std::map<size_t, size_t> rank_arbor;
        std::map<size_t, size_t> parent_arbor;
        // boost::disjoint_sets<AP_map, AP_map> dson_undirected_edges(rank, parent);
        boost::disjoint_sets<AP_map, AP_map> ds_arbor(
                boost::make_assoc_property_map(rank_arbor),
                boost::make_assoc_property_map(parent_arbor));

        for (size_t i = 0; i < n; ++i) {
            ds_contract.make_set(i);
            ds_arbor.make_set(i);
            parent_label_map_[i] = i;
        }
        Contract2(ds_contract);
        ExpandWithEdgeHandling();
    }

}
