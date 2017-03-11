#pragma once

#include <utils/fastq_reader.hpp>

struct QualityStatistics {
	double aver_read_qual;
	vector<double> aver_nucls_qual;

	void ShortPrint(ostream &out) {
		out << "Average read quality\t\t" << aver_read_qual << endl;
		out << "Average nucleotide quality:" << endl;
		for(auto it = aver_nucls_qual.begin(); it != aver_nucls_qual.end(); it++)
			out << *it << " ";
		out << endl;
	}

	void OutputNucleotideQuality(ostream &out) {
		for(auto it = aver_nucls_qual.begin(); it != aver_nucls_qual.end(); it++)
			out << *it << endl;
	}

	QualityStatistics() :
		aver_read_qual() { }
};

class QualityStatisticsCalculator {
	QualityStatistics stats_;
	size_t phred_offset_;
	size_t read_number_;
	size_t ideal_rl_;

	double AverageQuality(string quality) {
		size_t sum_qual = 0;
		for(size_t i = 0; i < quality.size(); i++)
			sum_qual += quality[i];
		return static_cast<double>(sum_qual) / static_cast<double>(quality.size());
	}

	void Initialize() {
		for(size_t i = 0; i < ideal_rl_; i++)
			stats_.aver_nucls_qual.push_back(0);
	}

public:
	QualityStatisticsCalculator(size_t phred_offset, size_t read_number, size_t ideal_rl) :
		phred_offset_(phred_offset),
		read_number_(read_number),
		ideal_rl_(ideal_rl) {
		Initialize();
	}

	void AddNuclsQuality(FastqRead read) {
		for(size_t i = 0; i < ideal_rl_; i++) {
			stats_.aver_nucls_qual[i] += double(size_t(read.quality[i])) / double(read_number_);
		}
	}

	void AddReadQuality(FastqRead read) {
		stats_.aver_read_qual += AverageQuality(read.quality) / static_cast<double>(read_number_);
	}

	QualityStatistics Stats() { return stats_; }

        size_t PhredOffset() const { return phred_offset_; }
};

class PairedReadQialityStatsCalculator {
	vector<PairedFastqRead> &reads_;
	QualityStatisticsCalculator calc_;

	size_t IdealReadLength() {
		size_t read_length = reads_[0].left_read.seq.size();
		for(size_t i = 0; i < reads_.size(); i++) {
			assert(reads_[i].left_read.seq.size() == read_length);
			assert(reads_[i].right_read.seq.size() == read_length);
		}
		//cout << "Ideal read length - " << read_length << endl;
		return read_length;
	}

public:
	PairedReadQialityStatsCalculator(vector<PairedFastqRead> &reads, size_t phred_offset) :
		reads_(reads),
		calc_(phred_offset, 2 * reads.size(), IdealReadLength()) {
	}

	void Calculate() {
		for(auto it = reads_.begin(); it != reads_.end(); it++) {
			calc_.AddReadQuality(it->left_read);
			calc_.AddReadQuality(it->right_read);
			calc_.AddNuclsQuality(it->left_read);
			calc_.AddNuclsQuality(it->right_read);
		}
	}

	QualityStatistics Stats() { return calc_.Stats(); }
};

class MergedReadQualityStatsCalculator {
	vector<FastqRead> &reads_;
	QualityStatisticsCalculator calc_;

	size_t ReadLength() {
		size_t rl(size_t(-1));
		for(auto it = reads_.begin(); it != reads_.end(); it++)
			rl = min<size_t>(it->seq.size(), rl);
		//cout << "Ideal read length - " << rl << endl;
		return rl;
	}

public:
	MergedReadQualityStatsCalculator(vector<FastqRead> &reads, size_t phred_offset) :
		reads_(reads),
		calc_(phred_offset, reads.size(), ReadLength()) { }

	void Calculate() {
		for(auto it = reads_.begin(); it != reads_.end(); it++) {
			calc_.AddNuclsQuality(*it);
			calc_.AddReadQuality(*it);
		}
	}

	QualityStatistics Stats() { return calc_.Stats(); }
};
