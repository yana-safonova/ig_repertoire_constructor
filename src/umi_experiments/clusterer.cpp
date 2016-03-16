#include "clusterer.hpp"

namespace clusterer {

    const ClusteringMode ClusteringMode::hamming = ClusteringMode(
            [](const seqan::Dna5String &first, const seqan::Dna5String &second) {
                return static_cast<size_t>(-half_sw_banded(first, second, 0, -1, -1, [](int) -> int { return 0; }, 0));
            }, 30);

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
    void GraphUmiPairsIterator::advance() {
        VERIFY_MSG(vertex_ < graph_->N(), "Trying to increment past-the-end iterator.");
        advances_ ++;

        ++ current_edge_;
        while (vertex_ < graph_->N() && current_edge_ == graph_->VertexEdges(vertex_).end()) {
            vertex_ ++;
            current_edge_ = graph_->VertexEdges(vertex_).begin();
            VERIFY_MSG(vertex_ < graph_->N() || advances_ == 2 * graph_->NZ(), "Expected to iterate over " << 2 * graph_->NZ() << " edges, but got " << advances_);
        }
    }

    GraphUmiPairsIterator GraphUmiPairsIterator::operator++() {
        advance();
        return *this;
    }

    GraphUmiPairsIterator GraphUmiPairsIterator::operator++(int) {
        const GraphUmiPairsIterator itr(*this);
        advance();
        return itr;
    }

    bool GraphUmiPairsIterator::operator==(GraphUmiPairsIterator other) const {
        VERIFY_MSG(graph_ == other.graph_, "Comparing iterators over different graph edges.");
        return vertex_ == other.vertex_ && current_edge_ == other.current_edge_;
    }

    bool GraphUmiPairsIterator::operator!=(GraphUmiPairsIterator other) const {
        return !(*this == other);
    }

    std::pair<size_t, size_t> GraphUmiPairsIterator::operator*() const {
        size_t opposite = *current_edge_;
        VERIFY_MSG(opposite < graph_->N(), "Bad vertex: " << opposite);
        return std::pair<size_t, size_t>(vertex_, opposite);
    }

    GraphUmiPairsIterable::GraphUmiPairsIterable(const SparseGraphPtr& graph) : graph_(graph), first_connected_vertex_(0) {
        while (first_connected_vertex_ < graph_->N() && graph_->Degree(first_connected_vertex_) == 0) {
            first_connected_vertex_ ++;
        }
    }

    GraphUmiPairsIterator GraphUmiPairsIterable::begin() const {
        return clusterer::GraphUmiPairsIterator(graph_, first_connected_vertex_, graph_->VertexEdges(first_connected_vertex_).begin());
    }

    GraphUmiPairsIterator GraphUmiPairsIterable::end() const {
        return clusterer::GraphUmiPairsIterator(graph_, graph_->N(), graph_->VertexEdges(graph_->N()).begin());
    }


    void FullGraphUmiPairsIterator::advance() {
        VERIFY_MSG(current_first_ < last_, "Trying to increment past-the-end iterator.");
        current_second_ ++;
        if ( current_second_ == current_first_ ) {
            current_first_ ++;
            current_second_ = 0;
        }
    }

    FullGraphUmiPairsIterator FullGraphUmiPairsIterator::operator++() {
        advance();
        return *this;
    }

    FullGraphUmiPairsIterator FullGraphUmiPairsIterator::operator++(int) {
        const FullGraphUmiPairsIterator itr(*this);
        advance();
        return itr;
    }

    bool FullGraphUmiPairsIterator::operator==(FullGraphUmiPairsIterator other) const {
        return current_first_ == other.current_first_ && current_second_ == other.current_second_;
    }

    bool FullGraphUmiPairsIterator::operator!=(FullGraphUmiPairsIterator other) const {
        return !(*this == other);
    }

    std::pair<size_t, size_t> FullGraphUmiPairsIterator::operator*() const {
        VERIFY_MSG(current_first_ < last_ && current_second_ < current_first_, "Dereferencing invalid iterator: " << current_first_ << ", " << current_second_);
        return std::pair<size_t, size_t>(current_first_, current_second_);
    }

    FullGraphUmiPairsIterator FullGraphUmiPairsIterable::begin() const {
        return clusterer::FullGraphUmiPairsIterator(1, 0, count_);
    }

    FullGraphUmiPairsIterator FullGraphUmiPairsIterable::end() const {
        return clusterer::FullGraphUmiPairsIterator(count_, 0, count_);
    }
}
