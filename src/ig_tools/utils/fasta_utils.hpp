#pragma once

struct fasta_read {
    std::string name;
    std::string seq;

    fasta_read() :
        name(),
        seq()
        { }

    fasta_read(std::string new_name, std::string new_seq) :
        name(new_name),
        seq(new_seq)
        { }

    void print(ostream &out) {
        out << "Name:\t" << name << std::endl;
        out << "Seq:\t" << seq << std::endl;
    }

    bool is_empty() {
        return name == "" || seq == "" ;
    }
};

class SingleFastaReader {
    ifstream src_;

public:
    SingleFastaReader(std::string fname) :
        src_(fname.c_str()) {
        VERIFY(!src_.fail());
    }
    fasta_read NextRead() {
        std::string name;
        std::string seq;
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

    std::vector<fasta_read> ReadFile() {
        std::vector<fasta_read> reads;
        while(!src_.eof()) {
            std::string name;
            std::string seq;
            getline(src_, name);
            getline(src_, seq);

            if(name != "" && seq != "")
                reads.push_back(fasta_read(name, seq));
        }
        return reads;
    }
};

class FastaWriter {
    std::ofstream out_;
public:
    FastaWriter(std::string fname) :
        out_(fname.c_str()) { }

    void Write(std::vector<fasta_read> reads) {
        for(size_t i = 0; i < reads.size(); i++) {
            out_ << reads[i].name << std::endl << reads[i].seq << std::endl <<
                    std::endl ;
        }
    }
    void Write(fasta_read read) {
        out_ << read.name << std::endl << read.seq << std::endl;
    }
};


