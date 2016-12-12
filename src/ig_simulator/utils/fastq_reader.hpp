#pragma once

#include "include_me.hpp"

struct fastq_read {
	string name;
	string seq;
	string quality;

	fastq_read() :
		name(),
		seq(),
		quality() { }

	fastq_read(string new_name, string new_seq, string new_quality) :
		name(new_name),
		seq(new_seq),
		quality(new_quality) { }

	void print(ostream &out) {
		out << "Name:\t" << name << endl;
		out << "Seq:\t" << seq << endl;
		out << "Qual:\t" << quality << endl;
	}

	bool is_empty() {
		return name == "" || seq == "" || quality == "";
	}
};

struct paired_fastq_read {
	fastq_read left_read;
	fastq_read right_read;

	paired_fastq_read() :
		left_read(),
		right_read() { }

	paired_fastq_read(fastq_read new_left, fastq_read new_right) :
		left_read(new_left),
		right_read(new_right) { }
};

class SingleFastqReader {
	ifstream src_;

public:
	SingleFastqReader(string fname) :
		src_(fname.c_str()) {
		assert(!src_.fail());
	}

	vector<fastq_read> ReadFile() {
		vector<fastq_read> reads;
		while(!src_.eof()) {
			string name;
			string seq;
			string tmp;
			string qual;
			getline(src_, name);
			getline(src_, seq);
			getline(src_, tmp);
			getline(src_, qual);

			assert(seq.size() == qual.size());

			if(name != "" && seq != "" && qual != "")
				reads.push_back(fastq_read(name, seq, qual));
		}
		return reads;
	}

    bool eof() {
        return src_.eof();
    }

    void close() {
        src_.close();
    }
    void reset(){
        src_.seekg(0, ios::beg);
    }
};

class PairedFastqReader {
	SingleFastqReader left_;
	SingleFastqReader right_;

public:
	PairedFastqReader(char *left_fname, char *right_fname) :
		left_(string(left_fname)),
		right_(string(right_fname)) {
	}

	vector<paired_fastq_read> Read() {
		vector<fastq_read> left_reads = left_.ReadFile();
		vector<fastq_read> right_reads = right_.ReadFile();
		assert(left_reads.size() == right_reads.size());
		vector<paired_fastq_read> paired_reads;
		for(size_t i = 0; i < left_reads.size(); i++)
			paired_reads.push_back(paired_fastq_read(left_reads[i],
					right_reads[i]));
		return paired_reads;
	}
};

class FastqWriter {
	ofstream out_;
public:
	FastqWriter(string fname) :
		out_(fname.c_str()) { }

	void Write(vector<fastq_read> reads) {
		for(size_t i = 0; i < reads.size(); i++) {
			out_ << reads[i].name << endl << reads[i].seq << endl <<
					"+" << endl << reads[i].quality << endl;
		}
	}
};
