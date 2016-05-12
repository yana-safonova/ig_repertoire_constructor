#pragma once

#include <sstream>
#include <vector>
#include <unordered_map>

namespace algorithms {
    struct Match {
        int subject_pos;
        int read_pos;
        size_t length;

        static int overlap(const Match &a,
                           const Match &b) {
            return std::max<int>(std::max<int>(a.length - (b.subject_pos - a.subject_pos),
                                               a.length - (b.read_pos - a.read_pos)), 0);
        }

        static bool less_subject_pos(const Match &a, const Match &b) {
            return a.subject_pos < b.subject_pos;
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
            return first().read_pos - first().subject_pos;
        }

        int right_shift() const {
            return last().read_pos - last().subject_pos;
        }

        int global_gap() const {
            return left_shift() - right_shift();
        }

        int read_segment_size() const {
            return last().read_pos - first().read_pos + last().length;
        }

        int subject_segment_size() const {
            return last().subject_pos - first().subject_pos + last().length;
        }

        std::string visualize_matches(int needle_length, int read_length) const;

        bool check_overlaps() const;
    };
}