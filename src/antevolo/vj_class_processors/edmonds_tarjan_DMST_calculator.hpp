#pragma once

#include <queue>
#include <stack>
#include <map>
#include <boost/property_map/property_map.hpp>
#include <boost/pending/disjoint_sets.hpp>

namespace antevolo {

    class EdmondsTarjanDMSTCalculator {
    public:
        struct WeightedEdge {
            size_t src_;
            size_t dst_;
            double weight_;

            WeightedEdge() :
                    src_(size_t(-1)),
                    dst_(size_t(-1)),
                    weight_(400) { }

            WeightedEdge(size_t src, size_t dst, double weight) :
                    src_(src),
                    dst_(dst),
                    weight_(weight) { }

            WeightedEdge(const WeightedEdge &oth) :
                    src_(oth.src_),
                    dst_(oth.dst_),
                    weight_(oth.weight_) { }

            WeightedEdge &operator=(const WeightedEdge& oth) {
                src_ = oth.src_;
                dst_ = oth.dst_;
                weight_ = oth.weight_;
                return *this;
            }
            bool operator== (const WeightedEdge& oth) const {
                return (src_ == oth.src_ && dst_ == oth.dst_ && weight_ == oth.weight_);
            }
            bool operator!= (const WeightedEdge& oth) const {
                return !(*this == oth);
            }
            bool operator> (const WeightedEdge& oth) const {
                return weight_ > oth.weight_;
            }
        };
        //struct Edge;
    private:


        typedef std::priority_queue<WeightedEdge,
                std::vector<WeightedEdge>,
                std::greater<WeightedEdge>> P_queue;
        typedef boost::associative_property_map<std::map<size_t, size_t>> AP_map;

        size_t n;
        size_t vertices_num_;
        size_t root_;
        const std::vector<WeightedEdge> &edges_;
        std::vector<WeightedEdge> in_;
        std::vector<double> const_add_;
        std::vector<size_t> parent_;
        std::vector<size_t> phase_;
        std::vector<std::vector<size_t>> children_;
        std::map<std::pair<size_t, size_t>, double> edge_initial_weight_map_;
        std::map<size_t, size_t> parent_label_map_;
        std::vector<P_queue> Ps_;
        std::stack<size_t> R_;

        void InitVertex();

        void Initialize();

        //void Contract(boost::disjoint_sets<AP_map, AP_map> &ds_contract,
        //              boost::disjoint_sets<AP_map, AP_map> &ds_arbor);

        void Contract2(boost::disjoint_sets<AP_map, AP_map> &ds_contract);

        //void Dismantle(size_t u);

        //void Expand(size_t r);

        void InEdgeHandle(size_t u);

        void ExpandWithEdgeHandling();


    public:
        /*
        struct Edge {
            size_t src_;
            size_t dst_;
            //double weight_;

            Edge() :
                    src_(size_t(-1)),
                    dst_(size_t(-1)) { }

            Edge(size_t src, size_t dst) :
                    src_(src),
                    dst_(dst) { }

            Edge(const Edge &oth) :
                    src_(oth.src_),
                    dst_(oth.dst_) { }

            Edge &operator=(Edge oth) {
                src_ = oth.src_;
                dst_ = oth.dst_;
                return *this;
            }

        };
        */



        EdmondsTarjanDMSTCalculator(size_t n, std::vector<WeightedEdge> &edges, size_t root) :
                n(n),
                vertices_num_(0),
                root_(root),
                edges_(edges) { }

        void EmpondsTarjan();

        std::vector<WeightedEdge> GetParentEdges() {
            std::vector<WeightedEdge> res;
            for (size_t i = 0; i < in_.size(); ++i) {
                WeightedEdge initial_e(in_[i]);
                initial_e.weight_ = edge_initial_weight_map_[std::make_pair(in_[i].src_, in_[i].dst_)];
                if (i != root_) {
                    res.push_back(initial_e);
                }
            }
            return res;
        }
        /*
        double GetEdgeWeight(const WeightedEdge& e)
        {
            return edge_weight_map_[e];
        }
        */

    };


}