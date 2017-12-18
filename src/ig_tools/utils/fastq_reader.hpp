#pragma once

#include <fstream>
#include <string>
#include <vector>
#include <verify.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

struct FastqRead {
	std::string name;
	std::string seq;
	std::string quality;

	FastqRead() :
		name(),
		seq(),
		quality() { }

	FastqRead(std::string new_name, std::string new_seq, std::string new_quality) :
		name(new_name),
		seq(new_seq),
		quality(new_quality) { }

	void print(std::ostream &out) {
		out << "Name:\t" << name << std::endl;
		out << "Seq:\t" << seq << std::endl;
		out << "Qual:\t" << quality << std::endl;
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
	std::ifstream file;
    std::istream src_;

    static bool has_suffix(const std::string &str, const std::string &suffix) {
        return str.size() >= suffix.size() &&
            str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    }

public:
	SingleFastqReader(std::string fname) : file(fname.c_str(), std::ios_base::in | std::ios_base::binary), src_(&fb) {
        if (has_suffix(fname, ".gz")) {
            fb.push(boost::iostreams::gzip_decompressor());
        } else {
            // Do nothing
        }

        fb.push(file);

		VERIFY(!src_.fail());
	}

	std::vector<FastqRead> ReadFile() {
		std::vector<FastqRead> reads;
		while(!src_.eof()) {
			std::string name;
			std::string seq;
			std::string tmp;
			std::string qual;
			getline(src_, name);
			getline(src_, seq);
			getline(src_, tmp);
			getline(src_, qual);

			VERIFY(seq.size() == qual.size());

			if(name != "" && seq != "" && qual != "")
				reads.push_back(FastqRead(name, seq, qual));
		}
		return reads;
	}

    bool eof() {
        return src_.eof();
    }

    void reset(){
        src_.seekg(0, std::ios::beg);
    }
};

class PairedFastqReader {
	SingleFastqReader left_;
	SingleFastqReader right_;

public:
	PairedFastqReader(char *left_fname, char *right_fname) :
		left_(std::string(left_fname)),
		right_(std::string(right_fname)) {
	}

	std::vector<PairedFastqRead> Read() {
		std::vector<FastqRead> left_reads = left_.ReadFile();
		std::vector<FastqRead> right_reads = right_.ReadFile();
		VERIFY(left_reads.size() == right_reads.size());
		std::vector<PairedFastqRead> paired_reads;
		for(size_t i = 0; i < left_reads.size(); i++)
			paired_reads.push_back(PairedFastqRead(left_reads[i],
					right_reads[i]));
		return paired_reads;
	}
};

class FastqWriter {
    std::ofstream out_;
public:
	FastqWriter(std::string fname) :
		out_(fname.c_str()) { }

	void Write(std::vector<FastqRead> reads) {
		for(size_t i = 0; i < reads.size(); i++) {
			out_ << reads[i].name << std::endl << reads[i].seq << std::endl <<
					"+" << std::endl << reads[i].quality << std::endl;
		}
	}
};
