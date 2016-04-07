#pragma once

#include <sstream>
#include <vector>
#include <unordered_map>

namespace algorithms {
    struct Match {
        int needle_pos;
        int read_pos;
        size_t length;

        static int overlap(const Match &a,
                           const Match &b) {
            return std::max<int>(std::max<int>(a.length - (b.needle_pos - a.needle_pos),
                                               a.length - (b.read_pos - a.read_pos)), 0);
        }

        static bool less_needle_pos(const Match &a, const Match &b) {
            return a.needle_pos < b.needle_pos;
        }

        static bool less_read_pos(const Match &a, const Match &b) {
            return a.read_pos < b.read_pos;
        }
    };

    class AlignmentPath : public std::vector<Match> {
        using std::vector<Match>::vector;
    public:
        int kplus_length() const;

        const Match &first() const {
            return (*this)[0];
        }

        const Match &last() const {
            return (*this)[this->size() - 1];
        }

        int left_shift() const {
            return first().read_pos - first().needle_pos;
        }

        int right_shift() const {
            return last().read_pos - last().needle_pos;
        }

        int global_gap() const {
            return left_shift() - right_shift();
        }

        int read_segment_size() const {
            return last().read_pos - first().read_pos + last().length;
        }

        int needle_segment_size() const {
            return last().needle_pos - first().needle_pos + last().length;
        }

        std::string visualize_matches(int needle_length, int read_length) const;

        bool check_overlaps() const;
    };
}