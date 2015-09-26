#pragma once

#include "../utils/include_me.hpp"

bool IsFastaName(string str) {
	if(str.size() == 0)
		return false;
	return str[0] == '>';
}

struct fasta_read {
	string name;
	string seq;

	fasta_read(string new_name, string new_seq) :
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
	ifstream &file_handler_;
public:
	PairedFastaReadConstructor(ifstream &file_handler) :
		file_handler_(file_handler) { }

	vector<paired_fasta_read> Read() {
		vector<paired_fasta_read> reads;
		string tmp_str;
		string name1, seq1, name2, seq2;
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
				assert(name1 != "");
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
		string tmp_str;
		string name, seq;
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
	FastaReader(string fname) :
		file_handler_(fname.c_str()),
		reader_(file_handler_) {
        if(file_handler_.fail()) {
            cout << "Input FASTA file " << fname << " was not found" << endl;
		    assert(!file_handler_.fail());
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
	ofstream out_;
public:
	FastaWriter(string fname) :
		out_(fname.c_str()) {
		assert(!out_.fail());
	}

	void Write(fasta_read read) {
		out_ << ">" << read.name << endl << read.seq << endl;
	}
};
