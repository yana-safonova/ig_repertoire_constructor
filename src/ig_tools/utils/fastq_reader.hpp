#pragma once

#include "include_me.hpp"
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

struct FastqRead {
	string name;
	string seq;
	string quality;

	FastqRead() :
		name(),
		seq(),
		quality() { }

	FastqRead(string new_name, string new_seq, string new_quality) :
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

struct PairedFastqRead {
	FastqRead left_read;
	FastqRead right_read;

	PairedFastqRead() :
		left_read(),
		right_read() { }

	PairedFastqRead(FastqRead new_left, FastqRead new_right) :
		left_read(new_left),
		right_read(new_right) { }
};

class SingleFastqReader {
    boost::iostreams::filtering_streambuf<boost::iostreams::input> fb;
    ifstream file;
    std::istream src_;

    static bool has_suffix(const std::string &str, const std::string &suffix) {
        return str.size() >= suffix.size() &&
            str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    }

public:
	SingleFastqReader(string fname) : file(fname.c_str(), ios_base::in | ios_base::binary), src_(&fb) {
        if (has_suffix(fname, ".gz")) {
            fb.push(boost::iostreams::gzip_decompressor());
        } else {
            // Do nothing
        }

        fb.push(file);

		assert(!src_.fail());
	}

	vector<FastqRead> ReadFile() {
		vector<FastqRead> reads;
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
				reads.push_back(FastqRead(name, seq, qual));
		}
		return reads;
	}

    bool eof() {
        return src_.eof();
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

	vector<PairedFastqRead> Read() {
		vector<FastqRead> left_reads = left_.ReadFile();
		vector<FastqRead> right_reads = right_.ReadFile();
		assert(left_reads.size() == right_reads.size());
		vector<PairedFastqRead> paired_reads;
		for(size_t i = 0; i < left_reads.size(); i++)
			paired_reads.push_back(PairedFastqRead(left_reads[i],
					right_reads[i]));
		return paired_reads;
	}
};

class FastqWriter {
	ofstream out_;
public:
	FastqWriter(string fname) :
		out_(fname.c_str()) { }

	void Write(vector<FastqRead> reads) {
		for(size_t i = 0; i < reads.size(); i++) {
			out_ << reads[i].name << endl << reads[i].seq << endl <<
					"+" << endl << reads[i].quality << endl;
		}
	}
};
