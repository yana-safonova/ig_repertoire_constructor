#include "clusterer.hpp"

namespace clusterer {

    ReadDist ClusteringMode::bounded_hamming_dist(size_t limit) {
        return [limit](const seqan::Dna5String& first, const seqan::Dna5String& second) {
            size_t diff = 0;
            size_t first_len = length(first);
            size_t second_len = length(second);
            size_t min_len = std::min(first_len, second_len);
            for (size_t i = 0; i < min_len; i ++) {
                if (first[i] != second[i]) {
                    diff ++;
                    if (diff > limit) {
                        return limit + 1;
                    }
                }
            }
            return diff + abs_diff(first_len, second_len);
        };
    }

    ReadDist ClusteringMode::bounded_edit_dist(size_t limit, size_t max_indels, bool binary) {
        return [limit, max_indels, binary](const seqan::Dna5String& first, const seqan::Dna5String& second) {
            const size_t INF = std::numeric_limits<size_t>::max() / 2;

            const size_t indel_cost = 1;
            const size_t mismatch_cost = 1;

            const size_t first_len = length(first);
            const size_t second_len = length(second);

            if (abs_diff(first_len, second_len) > std::min(limit, max_indels)) {
                return limit + 1;
            }

            std::vector<size_t> dp1(2 * max_indels + 1, INF);
            std::vector<size_t> dp2(2 * max_indels + 1);
            for (size_t j = 0; j <= max_indels && j <= second_len; j ++) {
                dp1[max_indels + j] = j;
            }

            for (size_t i = 0; i < first_len; i ++) {
                std::vector<size_t>& dp_cur = (i & 1) ? dp1 : dp2;
                std::vector<size_t>& dp_prev = (i & 1) ? dp2 : dp1;
                std::fill(dp_cur.begin(), dp_cur.end(), INF);

                for (size_t index = 0; index <= 2 * max_indels; index ++) {
                    if ( i + index < max_indels) continue;
                    if (i + index <= second_len + max_indels) {
                        if (index > 0) {
                            // i + 1 + index - 1 - max_indels
                            dp_cur[index - 1] = std::min(dp_cur[index - 1], dp_prev[index] + indel_cost);

                            // i + 1 + index (- 1) - max_indels
                            if (i + 1 + index <= second_len + max_indels) {
                                dp_cur[index] = std::min(dp_cur[index], dp_cur[index - 1] + indel_cost);
                            }
                        }

                        // i (+ 1) + index - max_indels
                        if (i + 1 + index <= second_len + max_indels) {
                            dp_cur[index] = std::min(dp_cur[index], dp_prev[index] + (first[i] == second[index - max_indels + i] ? 0 : mismatch_cost));
                        }
                    }
                }
//                INFO( i );
//                for (size_t index = 0; index <= 2 * max_indels; index ++) {
//                    std::cout << dp_cur[index] << " ";
//                }
//                std::cout << std::endl;

                bool all_too_large = true;
                for (size_t index = 0; index <= 2 * max_indels; index ++) {
                    size_t first_left = first_len - i - 1;
                    if (i + 1 + index > second_len + max_indels) break;
                    size_t second_left = second_len - i - 1 - index + max_indels;
                    size_t diff = abs_diff(first_left, second_left);
                    size_t lower = dp_cur[index] + diff;
                    if (lower <= limit) {
                        all_too_large = false;
                    }
                    size_t upper = dp_cur[index] + first_left + second_left;
                    if (binary && upper <= limit) {
//                        INFO("Pessimistic is " << upper << " for index " << index);
                        return upper;
                    }
                }
                if (all_too_large) {
//                    INFO("Optimistic is too high");
                    return limit + 1;
                }
            }
            std::vector<size_t>& dp = (first_len & 1) ? dp2 : dp1;
            return std::min(dp[max_indels + second_len - first_len], limit + 1);
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
