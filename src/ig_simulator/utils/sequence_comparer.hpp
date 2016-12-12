#pragma once

#include "include_me.hpp"
#include "sequence_tools.hpp"

struct SeqComparisonResults {
    bool match;
    bool second_reverse;
    pair<size_t, size_t> s1_overlap_pos;
    pair<size_t, size_t> s2_overlap_pos;

    SeqComparisonResults():
        match(false),
        second_reverse(false),
        s1_overlap_pos(make_pair(0, 0)),
        s2_overlap_pos(make_pair(0, 0)) { }

    SeqComparisonResults(pair<size_t, size_t> new_s1_pos, pair<size_t, size_t> new_s2_pos):
        match(true),
        second_reverse(false),
        s1_overlap_pos(new_s1_pos),
        s2_overlap_pos(new_s2_pos) { }

};

class SequenceComparer {
	static const size_t min_overlap = 200;

	SeqComparisonResults TwoSequencesMatch(string s1, string s2) {
		pair<size_t, size_t> start_pos(0, s2.size() - 1);
		for(size_t i = 0; i < s1.size() + s2.size() - 1; i++) {
			size_t overlap_size = min<size_t>(s1.size() - start_pos.first,
					s2.size() - start_pos.second);
			if(overlap_size >= min_overlap) {
				bool sequences_match = true;
				for(size_t j = 0; j < overlap_size; j++)
					if(s1[start_pos.first + j] != s2[start_pos.second + j]) {
						sequences_match = false;
						break;
					}
				if(sequences_match) {
                    return SeqComparisonResults(
                        make_pair(start_pos.first, start_pos.first + overlap_size - 1),
                                    make_pair(start_pos.second, start_pos.second + overlap_size - 1));
				}
			}
			if(start_pos.second > 0 && start_pos.first == 0) {
				start_pos.second--;
			} else if(start_pos.second == 0)
				start_pos.first++;
		}
		return SeqComparisonResults();
	}

public:
	SeqComparisonResults SequencesMatch(string s1, string s2) {
		string rc_s2 = reverse_complementary(s2);
		auto direct_res = TwoSequencesMatch(s1, s2);
        if(direct_res.match) {
            direct_res.second_reverse = false;
            return direct_res;
        }

        auto reverse_res = TwoSequencesMatch(s1, rc_s2);
        if(reverse_res.match) {
            reverse_res.second_reverse = true;
            return reverse_res;
        }

        return SeqComparisonResults();
	}
};

class MismatchSequenceComparer {
    static const size_t min_overlap = 200;
    size_t max_dist_;

    size_t GetBestHammingDistance(string s1, string s2) {
		pair<size_t, size_t> start_pos(0, s2.size() - 1);
		for(size_t i = 0; i < s1.size() + s2.size() - 1; i++) {
			size_t overlap_size = min<size_t>(s1.size() - start_pos.first,
					s2.size() - start_pos.second);
			if(overlap_size >= min_overlap) {
                //string substr1 = s1.substr(start_pos.first, overlap_size);
                //string substr2 = s2.substr(start_pos.second, overlap_size);

                size_t cur_dist = 0;
                for(size_t i = 0; i < overlap_size; i++)
                    if(s1[start_pos.first + i] != s2[start_pos.second + i]) {
                        cur_dist++;
                        if(cur_dist > max_dist_)
                            continue;
                    }
                        
                //size_t dist = HammingDistance(substr1, substr2);
                if(cur_dist <= max_dist_)
                    return cur_dist;
            }
			if(start_pos.second > 0 && start_pos.first == 0) {
				start_pos.second--;
			} else if(start_pos.second == 0)
				start_pos.first++;
		}
        return min_overlap;
    }

public:
    MismatchSequenceComparer(size_t max_dist) :
        max_dist_(max_dist) { }

    size_t CompareString(string s1, string s2) {
        string rc_s2 = reverse_complementary(s2);
        return min<size_t>(GetBestHammingDistance(s1, s2), GetBestHammingDistance(s1, rc_s2)); 
    }
};





















