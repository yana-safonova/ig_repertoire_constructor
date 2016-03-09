#include "utils.hpp"

size_t get_sw_dist(const seqan::Dna5String& first, const seqan::Dna5String& second) {
    std::vector<size_t> cur(length(first) + 1);
    std::vector<size_t> prev(length(first) + 1);
    size_t INF = std::numeric_limits<size_t>::max() / 2;
    std::fill(cur.begin() + 1, cur.end(), INF);
    for (size_t i = 0; i < length(second); i ++) {
        std::swap(cur, prev);
        std::fill(cur.begin(), cur.end(), INF);
        cur[0] = prev[0] + 1;
        for (size_t j = 1; j <= length(first); j ++) {
            cur[j] = std::min({
                                 cur[j - 1] + 1,
                                 prev[j] + 1,
                                 prev[j - 1] + (first[j - 1] == second[i] ? 0 : 1)
                         });
        }
    }
    return cur[length(first)];
}

void create_console_logger() {
    using namespace logging;
    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}
