#include <verify.hpp>
#include "block_alignment_primitives.hpp"

#include <boost/format.hpp>

namespace algorithms {
    size_t AlignmentPath::kplus_length() const {
        size_t result = 0;
        for (const auto &match : *this)
            result += match.length;
        return result;
    }

    std::string AlignmentPath::visualize_matches(int needle_length, int read_length) const {
        const std::vector <Match> &matches = *this; // TODO Fixit

        // Draw fancy alignment
        // (read/needle)
        std::stringstream ss;

        VERIFY(std::is_sorted(matches.cbegin(), matches.cend(), Match::less_subject_pos));
        VERIFY(std::is_sorted(matches.cbegin(), matches.cend(), Match::less_read_pos));

        ss << boost::format("{%d}") % std::max(matches[0].subject_pos - matches[0].read_pos, 0);
        ss << boost::format("(%d)") % std::min(matches[0].subject_pos, matches[0].read_pos);
        for (size_t i = 0; i < matches.size() - 1; ++i) {
            int read_gap = matches[i + 1].read_pos - matches[i].read_pos - int(matches[i].length);
            int needle_gap = matches[i + 1].subject_pos - matches[i].subject_pos - int(matches[i].length);
            size_t current_match_len = matches[i].length;

            if (needle_gap >= 0 || read_gap >= 0) {
                if (needle_gap != read_gap) {
                    ss << boost::format("%d(%d%+d)") % current_match_len % needle_gap % (read_gap - needle_gap);
                } else {
                    ss << boost::format("%2%(%1%)") % read_gap % current_match_len;
                }
            }
        }

        const auto &last_match = matches[matches.size() - 1];
        ss << boost::format("%1%") % last_match.length;
        ss << boost::format("(%d)") % std::min(needle_length - last_match.subject_pos - last_match.length,
                                         read_length - last_match.read_pos - last_match.length);
        ss << boost::format("{%d}") %
              std::max((needle_length - last_match.subject_pos) - (read_length - last_match.read_pos), 0);

        return ss.str();
    }

    bool AlignmentPath::check_overlaps() const {
        if (size() > 1) {
            for (size_t i = 0; i < size() - 1; ++i) {
                if (Match::overlap((*this)[i], (*this)[i + 1])) {
                    return false;
                }
            }
        }
        return true;
    }
}