#pragma once

#include <string>
#include <vector>
#include <fstream>

bool IsFastaName(std::string str) {
	if(str.size() == 0)
		return false;
	return str[0] == '>';
}

struct fasta_read {
	std::string name;
	std::string seq;

	fasta_read(std::string new_name, std::string new_seq) :
		name(new_name),
		seq(new_seq) { }

	fasta_read() :
		name(), seq() { }
};

struct paired_fasta_read {
	fasta_read left_read;
	fasta_read right_read;

	paired_fasta_read() :
		left_read(), right_read() { }

	paired_fasta_read(fasta_read new_left_read,
			fasta_read new_right_read) :
		left_read(new_left_read),
		right_read(new_right_read) { }
};

class PairedFastaReadConstructor {
	std::ifstream &file_handler_;
public:
	PairedFastaReadConstructor(std::ifstream &file_handler) :
		file_handler_(file_handler) { }

    std::vector<paired_fasta_read> Read() {
        std::vector<paired_fasta_read> reads;
		std::string tmp_str;
		std::string name1, seq1, name2, seq2;
		while(!file_handler_.eof()) {
			getline(file_handler_, tmp_str);
			if(IsFastaName(tmp_str)){
				if(name1 == "" && name2 == "" && seq1 == "" && seq2 == "")
					name1 = tmp_str;
				else
					if(name1 != "" && name2 != "" && seq1 != "" && seq2 != "") {
						reads.push_back(paired_fasta_read(fasta_read(name1, seq1), fasta_read(name2, seq2)));
						name1 = tmp_str; name2 = ""; seq1 = ""; seq2 = "";
					}
					else
						name2 = tmp_str;
			}
			else{
				VERIFY(name1 != "");
				if(name2 == "")
					seq1 = seq1 + tmp_str;
				else
					seq2 = seq2 + tmp_str;
			}
		}
		if(name1 != "" && name2 != "" && seq1 != "" && seq2 != "")
			reads.push_back(paired_fasta_read(fasta_read(name1, seq1), fasta_read(name2, seq2)));
		return reads;
	}
};

class FastaReadConstructor {
	ifstream &file_handler_;
public:
	FastaReadConstructor(ifstream &file_handler) :
		file_handler_(file_handler) { }

	vector<fasta_read> Read() {
		vector<fasta_read> reads;
		std::string tmp_str;
		std::string name, seq;
		while(!file_handler_.eof()) {
			getline(file_handler_, tmp_str);
			if(tmp_str != "") {
				if(IsFastaName(tmp_str)){
					if(name != "" && seq != "")
						reads.push_back(fasta_read(name, seq));
					name = tmp_str;
					seq = "";
				}
				else
						seq = seq + tmp_str;
			}
		}
		if(name != "" && seq != "")
			reads.push_back(fasta_read(name, seq));
		return reads;
	}
};

template<class ReadType, class ReadConstructor>
class FastaReader {
	ifstream file_handler_;
	ReadConstructor reader_;
public:
	FastaReader(std::string fname) :
		file_handler_(fname.c_str()),
		reader_(file_handler_) {
        if(file_handler_.fail()) {
            cout << "Input FASTA file " << fname << " was not found" << std::endl;
			VERIFY(!file_handler_.fail());
        }
	}

	vector<ReadType> Read() {
		return reader_.Read();
	}

	void Release() {
		file_handler_.close();
	}

	~FastaReader() {
		file_handler_.close();
	}
};

class FastaWriter {
	std::ofstream out_;
public:
	FastaWriter(std::string fname) :
		out_(fname.c_str()) {
		VERIFY(!out_.fail());
	}

	void Write(fasta_read read) {
		out_ << ">" << read.name << std::endl << read.seq << std::endl;
	}
};
