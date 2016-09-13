#include "clusterer.hpp"

namespace clusterer {

    ReadDist ClusteringMode::bounded_hamming_dist(size_t limit) {
        return [limit](const seqan::Dna5String& first, const seqan::Dna5String& second) {
            size_t diff = 0;
            size_t len = std::min(length(first), length(second));
            for (size_t i = 0; i < len; i ++) {
                if (first[i] != second[i]) {
                    diff ++;
                    if (diff > limit) {
                        return limit + 1;
                    }
                }
            }
            return diff;
        };
    }

    ReadDist ClusteringMode::bounded_edit_dist(size_t limit, size_t max_indels, bool binary) {
        return [limit, max_indels, binary](const seqan::Dna5String& first, const seqan::Dna5String& second) {
            const size_t INF = std::numeric_limits<size_t>::max() / 2;

            const size_t indel_cost = 1;
            const size_t mismatch_cost = 1;

            const size_t first_len = length(first);
            const size_t second_len = length(second);

            std::vector<size_t> dp_cur(2 * max_indels + 1, INF);
            std::vector<size_t> dp_prev(2 * max_indels + 1);
            for (size_t j = 0; j <= max_indels && j <= second_len; j ++) {
                dp_cur[max_indels + j] = j;
            }
            // lizard tail for first sequence: minimum of values edit_dist(prefix of first, second)
            size_t min_at_top = second_len <= max_indels ? dp_cur[max_indels + second_len] : INF;

            for (size_t i = 0; i < first_len; i ++) {
                std::swap(dp_cur, dp_prev);
                std::fill(dp_cur.begin(), dp_cur.end(), INF);

                for (size_t index = 0; index <= 2 * max_indels; index ++) {
                    if (index > 0) {
                        dp_cur[index - 1] = std::min(dp_cur[index - 1], dp_prev[index] + indel_cost);
                        dp_cur[index] = std::min(dp_cur[index], dp_cur[index - 1] + indel_cost);
                    }
                    if (index + i >= max_indels && index - max_indels + i < second_len) {
                        dp_cur[index] = std::min(dp_cur[index], dp_prev[index] + (first[i] == second[index - max_indels + i] ? 0 : mismatch_cost));
                    }
                }

                if (second_len + max_indels >= i + 1 && second_len - i - 1 + max_indels < 2 * max_indels + 1) {
                    min_at_top = std::min(min_at_top, dp_cur[second_len - i - 1 + max_indels]);
                    if (binary && min_at_top <= limit) return min_at_top;
                }
                if (*std::max_element(dp_cur.begin(), dp_cur.end()) > limit) return limit + 1;
            }
            size_t result = std::min(min_at_top, *std::min_element(dp_cur.begin(), dp_cur.end()));

            return std::min(result, limit + 1);
        };
    }

    ClusterDistChecker ClusteringMode::clusters_close_by_center(const ReadDist& check_read_dist, size_t limit) {
        return [check_read_dist, limit](const ClusterPtr<Read>& first, const ClusterPtr<Read>& second) {
            return check_read_dist(first->center, second->center) <= limit;
        };
    }

    ClusterDistChecker ClusteringMode::clusters_close_by_min(const ReadDist& read_dist, size_t limit) {
        return [read_dist, limit](const ClusterPtr<Read>& first, const ClusterPtr<Read>& second) {
            for (const auto& first_read : first->members) {
                for (const auto& second_read : second->members) {
                    if (read_dist(first_read.GetSequence(), second_read.GetSequence()) <= limit) {
                        return true;
                    }
                }
            }
            return false;
        };
    }


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
