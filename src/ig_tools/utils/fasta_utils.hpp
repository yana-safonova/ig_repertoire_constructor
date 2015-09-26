#pragma once

#include "include_me.hpp"

struct fasta_read {
    string name;
    string seq;

    fasta_read() :
        name(),
        seq()
        { }

    fasta_read(string new_name, string new_seq) :
        name(new_name),
        seq(new_seq)
        { }

    void print(ostream &out) {
        out << "Name:\t" << name << endl;
        out << "Seq:\t" << seq << endl;
    }

    bool is_empty() {
        return name == "" || seq == "" ;
    }
};

class SingleFastaReader {
    ifstream src_;

public:
    SingleFastaReader(string fname) :
        src_(fname.c_str()) {
        assert(!src_.fail());
    }
    fasta_read NextRead() {
        string name;
        string seq;
        getline(src_, name);
        getline(src_, seq);
//        cout << name<< " "<< seq.size() << " " << qual.size() << '\n';

        return fasta_read(name, seq) ;
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

    vector<fasta_read> ReadFile() {
        vector<fasta_read> reads;
        while(!src_.eof()) {
            string name;
            string seq;
            getline(src_, name);
            getline(src_, seq);

            if(name != "" && seq != "")
                reads.push_back(fasta_read(name, seq));
        }
        return reads;
    }
};

class FastaWriter {
    ofstream out_;
public:
    FastaWriter(string fname) :
        out_(fname.c_str()) { }

    void Write(vector<fasta_read> reads) {
        for(size_t i = 0; i < reads.size(); i++) {
            out_ << reads[i].name << endl << reads[i].seq << endl <<
                    endl ;
        }
    }
    void Write(fasta_read read) {
        out_ << read.name << endl << read.seq << endl;
    }
};


