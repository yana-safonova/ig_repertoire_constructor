#pragma once

#include "../utils/fastq_reader.hpp"
#include "../utils/include_me.hpp"
#include "../utils/sequence_tools.hpp"
#include "../utils/string_tools.hpp"

struct merger_setting {
	size_t min_overlap;
	double max_mismatch_rate;
	bool simulated_mode;

	merger_setting() :
		min_overlap(50),
		max_mismatch_rate(.1),
		simulated_mode(false) { }

	void print(ostream &out) {
		out << "Sequence merger settings:" << endl;
		out << "Min overlap size\t" << min_overlap << endl;
		out << "Max mismatch rate\t" << max_mismatch_rate << endl;
		if(simulated_mode)
			out << "Simulated mode is ON" << endl;
	}
};

class SequenceMerger {
	merger_setting setting_;

	string reverse_quality_string(string qual) {
		for(size_t i = 0; i < qual.size() / 2; i++) {
			char tmp = qual[i];
			qual[i] = qual[qual.size() - i - 1]; 
			qual[qual.size() - i - 1] = tmp;
		}
		return qual;
	}

	pair<size_t, size_t> FindBestOverlap(string seq1, string seq2) {
		pair<size_t, size_t> best_overlap(size_t(-1), size_t(-1));
		for(size_t i = 0; i < seq1.size(); i++) {
			size_t overlap_size = min<size_t>(seq2.size(), seq1.size() - i);
			if(overlap_size >= setting_.min_overlap) {
				string subseq1 = seq1.substr(i, overlap_size);
				string subseq2 = seq2.substr(0, overlap_size);
				size_t dist = HammingDistance(subseq1, subseq2);
				double mism_rate = double(dist) / overlap_size;
				if(mism_rate <= setting_.max_mismatch_rate) {
					if(best_overlap.first > dist) {
						best_overlap.first = dist;
						best_overlap.second = i;
					}
				}
			}
		}
		return best_overlap;
	}

	pair<string, string> MergeQualifiedSeq(paired_fastq_read &paired_read) {
		string rc_right = reverse_complementary(paired_read.right_read.seq);
		string left = paired_read.left_read.seq;
		pair<size_t, size_t> overlap = FindBestOverlap(left, rc_right);
		if(overlap == make_pair(size_t(-1), size_t(-1)))
			return make_pair(string(), string());

		size_t overlap_size = min<size_t>(left.size() - overlap.second,
				rc_right.size());
		string overlap_seq = left.substr(overlap.second, overlap_size);
		string qual1 = paired_read.left_read.quality;
		string qual2 = reverse_quality_string(paired_read.right_read.quality);
		string overlap_qual = qual1.substr(overlap.second,
				overlap_size);

		size_t mismatch_number = 0;
		if(!setting_.simulated_mode)
			for(size_t i = 0; i < overlap_size; i++)
				if(qual1[overlap.second + i] < qual2[i]) {
					overlap_seq[i] = rc_right[i];
					overlap_qual[i] = qual2[i];
				}

		string merged_seq = left.substr(0, overlap.second) + overlap_seq;
		string merged_qual = qual1.substr(0, overlap.second) + overlap_qual;
		if(overlap_size > rc_right.size()) {
			merged_seq = merged_seq + left.substr(overlap.second + overlap_size,
					left.size() - overlap.second - overlap_size);
			merged_qual = merged_qual + qual1.substr(overlap.second + overlap_size,
					left.size() - overlap.second - overlap_size);
		}
		else {
			merged_seq = merged_seq + rc_right.substr(overlap_size,
					rc_right.size() - overlap_size);
			merged_qual = merged_qual + qual2.substr(overlap_size,
					rc_right.size() - overlap_size );
		}
		assert(merged_seq.size() == merged_qual.size());
		return make_pair(merged_seq, merged_qual);
	}

	string MergeNames(size_t index, string name1, string name2) {
		if(name1[0] == '@')
			name1 = name1.substr(1, name1.size() - 1);
		stringstream ss;
		ss << "@" << index << "_merged_read_" << delete_spaces(name1);
		return ss.str();
	}

public:
	SequenceMerger(merger_setting setting) :
		setting_(setting) { }

	fastq_read Merge(size_t index, paired_fastq_read paired_read) {
		string merged_name = MergeNames(index, paired_read.left_read.name,
				paired_read.right_read.name);
		pair<string, string> merged_seq_qual = MergeQualifiedSeq(paired_read);
		return fastq_read(merged_name, merged_seq_qual.first,
				merged_seq_qual.second);
	}
};

class PairedReadsMerger {
	SequenceMerger seq_merger_;
public:
	PairedReadsMerger(merger_setting settings) :
		seq_merger_(settings) { }

	vector<fastq_read> Merge(vector<paired_fastq_read> reads) {
		vector<fastq_read> merged_reads;
		double processed_rate = 10;
		for(size_t i = 0; i < reads.size(); i++) {
			fastq_read merged_read = seq_merger_.Merge(i, reads[i]);
			if(!merged_read.is_empty())
				merged_reads.push_back(merged_read);
			if(double(i) / reads.size() * 100 > processed_rate) {
				cout << processed_rate << "% reads were processed" << endl;
				processed_rate += 10;
			}
		}
		cout << "100% reads were processed" << endl;
		return merged_reads;
	}
};
