#include "clusterer.hpp"

namespace clusterer {

    const ClusteringMode ClusteringMode::hamming = ClusteringMode(
            [](const seqan::Dna5String &first, const seqan::Dna5String &second) {
                return static_cast<size_t>(-half_sw_banded(first, second, 0, -1, -1, [](int) -> int { return 0; }, 0));
            }, 10);

    const ClusteringMode ClusteringMode::edit = ClusteringMode(
            [](const seqan::Dna5String &first, const seqan::Dna5String &second) {
                return static_cast<unsigned long>(get_sw_dist(first, second));
            }, 10);


    ReflexiveUmiPairsIterator ReflexiveUmiPairsIterator::operator++() {
        current_ ++;
        return *this;
    }

    ReflexiveUmiPairsIterator ReflexiveUmiPairsIterator::operator++(int) {
        const ReflexiveUmiPairsIterator itr(*this);
        current_ ++;
        return itr;
    }

    std::pair<size_t, size_t> ReflexiveUmiPairsIterator::operator*() const {
        return std::make_pair(current_, current_);
    }

    bool ReflexiveUmiPairsIterator::operator==(ReflexiveUmiPairsIterator other) const {
        VERIFY_MSG(other.last_ == last_, "Comparing different iterators.");
        return current_ == other.current_;
    }

    bool ReflexiveUmiPairsIterator::operator!=(ReflexiveUmiPairsIterator other) const {
        return !(*this == other);
    }

    ReflexiveUmiPairsIterator ReflexiveUmiPairsIterable::begin() const {
        return ReflexiveUmiPairsIterator(0, count_);
    }

    ReflexiveUmiPairsIterator ReflexiveUmiPairsIterable::end() const {
        return ReflexiveUmiPairsIterator(count_, count_);
    }

    // TODO: skip half of the edges
    GraphUmiPairsIterator GraphUmiPairsIterator::operator++() {
        VERIFY_MSG(vertex_ < graph_->N(), "Trying to increment past-the-end iterator.");
        advances_ ++;

        ++ *current_edge_;
        while (vertex_ < graph_->N() && *current_edge_ == graph_->VertexEdges(vertex_).end()) {
            vertex_ ++;
            // euw, ugly.. how to do this??
            current_edge_ = std::make_shared<SparseGraph::EdgesIterator>(graph_->VertexEdges(vertex_).begin());
            VERIFY_MSG(vertex_ < graph_->N() || advances_ == 2 * graph_->NZ(), "Expected to iterate over " << 2 * graph_->NZ() << " edges, but got " << advances_);
        }
        return *this;
    }

    GraphUmiPairsIterator GraphUmiPairsIterator::operator++(int) {
        VERIFY_MSG(vertex_ < graph_->N(), "Trying to increment past-the-end iterator.");
        advances_ ++;

        const GraphUmiPairsIterator itr(*this);
        ++ *current_edge_;
        while (vertex_ < graph_->N() && *current_edge_ == graph_->VertexEdges(vertex_).end()) {
            vertex_ ++;
            current_edge_ = std::make_shared<SparseGraph::EdgesIterator>(graph_->VertexEdges(vertex_).begin());
            VERIFY_MSG(vertex_ < graph_->N() || advances_ == 2 * graph_->NZ(), "Expected to iterate over " << 2 * graph_->NZ() << " edges, but got " << advances_);
        }
        return itr;
    }

    bool GraphUmiPairsIterator::operator==(GraphUmiPairsIterator other) const {
        VERIFY_MSG(graph_ == other.graph_, "Comparing iterators over different graph edges.");
        return vertex_ == other.vertex_ && *current_edge_ == *other.current_edge_;
    }

    bool GraphUmiPairsIterator::operator!=(GraphUmiPairsIterator other) const {
        return !(*this == other);
    }

    std::pair<size_t, size_t> GraphUmiPairsIterator::operator*() const {
        size_t opposite = **current_edge_;
        VERIFY_MSG(opposite < graph_->N(), "Bad vertex: " << opposite);
        return std::pair<size_t, size_t>(vertex_, opposite);
    }

    GraphUmiPairsIterator GraphUmiPairsIterable::begin() const {
        return clusterer::GraphUmiPairsIterator(graph_, 0, graph_->VertexEdges(0).begin());
    }

    GraphUmiPairsIterator GraphUmiPairsIterable::end() const {
        return clusterer::GraphUmiPairsIterator(graph_, graph_->N(), graph_->VertexEdges(graph_->N()).begin());
    }
}
