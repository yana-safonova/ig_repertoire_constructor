#pragma once

#include <fstream>
#include <string>
#include <vector>
#include "../fasta_reader.hpp"
#include "../string_tools.hpp"

struct SeqCluster {
    std::string seq;
    size_t size;
    size_t id;
};

class SeqClusterParser {
public:
    static SeqCluster Parse(fasta_read read) {
        SeqCluster cluster;
        cluster.seq = read.seq;
        std::vector<std::string> header_splits = split(read.name, "___");
        cluster.id = string_to_number<size_t>(header_splits[1]);
        cluster.size = string_to_number<size_t>(header_splits[3]);
        return cluster;
    }
};

//class SeqClusterization {
//    std::vector<SeqClass> clusters_;
//};

class ClustersFastaReader {
    FastaReader<fasta_read, FastaReadConstructor> fasta_reader_;
    std::string fname_;
public:
    ClustersFastaReader(std::string fname):
        fasta_reader_(fname),
        fname_(fname) { }

    std::vector<SeqCluster> Read() {
        std::vector<fasta_read> reads = fasta_reader_.Read();
        std::vector<SeqCluster> clusters;
        for(auto it = reads.begin(); it != reads.end(); it++) {
            clusters.push_back(SeqClusterParser::Parse(*it));
        }
        std::cout << clusters.size() << " cluster(s) were extracted from " << fname_ << std::endl;
        return clusters;
    }
};

struct stop_codon_pos {
        bool valid;
        size_t pos1;
        size_t pos2;
        size_t pos3;

        void SetValid(size_t seq_size) { valid = pos1 < seq_size && pos2 < seq_size && pos3 < seq_size; }

        bool operator==(const stop_codon_pos& obj) {
                return pos1 == obj.pos1 && pos2 == obj.pos2 && pos3 && obj.pos3;
        }
};


class StopCodonSearcher {
        static size_t StopCodonOccurs(std::string seq, size_t shift) {
                std::string stop_codon1 = "TAG";
                std::string stop_codon2 = "TAA";
                std::string stop_codon3 = "TGA";

                size_t start_pos = 3;
                size_t end_pos = seq.size() - 3;

                size_t start_aa = start_pos / 3 + (start_pos % 3 != 0) ? 1 : 0;
                size_t num_aa = (end_pos - shift) / 3;

                for(size_t i = start_aa; i < num_aa; i++) {
                        size_t pos = i * 3 + shift;
                        std::string aa = seq.substr(pos, 3);
                        if(aa == stop_codon1 || aa == stop_codon2 || aa == stop_codon3)
                                return pos;
                }
                return seq.size();
        }

public:
        static stop_codon_pos SequenceContainsStopCodon(std::string seq) {
                stop_codon_pos res;
                res.pos1 = StopCodonOccurs(seq, 0);
                res.pos2 = StopCodonOccurs(seq, 1);
                res.pos3 = StopCodonOccurs(seq, 2);
                res.SetValid(seq.size());
                return res;
        }
};


class AAVerificator {
    std::vector<size_t> incorrect_clusters_ids_;
    std::vector<size_t> incorrect_clusters_sizes_;
public:
    void Verify(std::vector<SeqCluster> &clusters) {
        for(auto it = clusters.begin(); it != clusters.end(); it++) {
            if(StopCodonSearcher::SequenceContainsStopCodon(it->seq).valid) {
                incorrect_clusters_ids_.push_back(it->id);
                incorrect_clusters_sizes_.push_back(it->size);
                //std::cout << it->id << "\t" << it->seq << std::endl;
            }
        }
        std::cout << incorrect_clusters_ids_.size() << " cluster contain stop codons" << std::endl;
    }

    void WriteIdsToFile(std::string output_fname) {
        std::ofstream out(output_fname.c_str());
        for(auto it = incorrect_clusters_ids_.begin(); it != incorrect_clusters_ids_.end(); it++)
            out << *it << std::endl;
    }

    void WriteSizesToFile(std::string output_fname) {
        std::ofstream out(output_fname.c_str());
        for(auto it = incorrect_clusters_sizes_.begin(); it != incorrect_clusters_sizes_.end(); it++)
            out << *it << std::endl;
    }
};
