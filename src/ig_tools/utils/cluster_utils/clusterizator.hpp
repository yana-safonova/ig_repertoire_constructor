#pragma once

#include "../include_me.hpp"
#include "../concurrent_dsu.hpp"
#include "../sequence_comparer.hpp"
#include "../fastq_reader.hpp"

template<class ReadType>
class ReadIndex {
protected:
	const vector<ReadType> &reads_;
public:
	ReadIndex(const vector<ReadType> &reads) :
		reads_(reads) {
	}

	virtual ~ReadIndex() {
	}

	virtual set<size_t> GetCandidatesFor(size_t index) = 0;
};

template<class ReadType>
class EasyReadIndex : public ReadIndex<ReadType> {
	vector<size_t> indices_;
public:
	EasyReadIndex(const vector<ReadType> &reads) :
		ReadIndex<ReadType>(reads) {
		for(size_t i = 0; i < ReadIndex<ReadType>::reads_.size(); i++)
			indices_.push_back(i);
	}

	set<size_t> GetCandidatesFor(size_t index) {
		return set<size_t>(indices_.begin() + index + 1, indices_.end());
	}
};

template<class ReadType>
class FLReadIndex : public ReadIndex<ReadType> {
	map<string, vector<size_t> > index_;
	const size_t k_;
public:
	FLReadIndex(const vector<ReadType> &reads, size_t k) :
		ReadIndex<ReadType>(reads), k_(k) {
		for(size_t i = 0; i < ReadIndex<ReadType>::reads_.size(); i++) {
			if(ReadIndex<ReadType>::reads_[i].seq.size() >= k) {
				string start_kmer = ReadIndex<ReadType>::reads_[i].seq.substr(0, k);
				string end_kmer = ReadIndex<ReadType>::reads_[i].seq.substr(ReadIndex<ReadType>::reads_[i].seq.size() - k, k);

				index_[start_kmer].push_back(i);
				index_[end_kmer].push_back(i);
				index_[reverse_complementary(start_kmer)].push_back(i);
				index_[reverse_complementary(end_kmer)].push_back(i);
			}
		}

        //for(auto it = index_.begin(); it != index_.end(); it++) {
        //    cout << it->first << ": ";
        //    for(auto it2 = it->second.begin(); it2 != it->second.end(); it2++) 
        //        cout << *it2 << " ";
        //    cout << endl;
        //}
        
	}

	set<size_t> GetCandidatesFor(size_t index) {
		set<size_t> result;
		for(size_t i = 0; i + k_ <= ReadIndex<ReadType>::reads_[index].seq.size(); i++) {
			string kmer = ReadIndex<ReadType>::reads_[index].seq.substr(i, k_);

            //cout << "Kmer : " << kmer << endl;

			auto it = index_.find(kmer);
			if(it != index_.end()) {
                for(auto it2 = it->second.begin(); it2 != it->second.end(); it2++)
                    result.insert(*it2);
				//result.insert(result.end(), it->second.begin(), it->second.end());
			}

			string rc_kmer = reverse_complementary(kmer);
			it = index_.find(rc_kmer);
			if(it != index_.end()) {
                for(auto it2 = it->second.begin(); it2 != it->second.end(); it2++)
                    result.insert(*it2);
				//result.insert(result.end(), it->second.begin(), it->second.end());
			}
		}

        //assert(false);
		return result;
	}
};

class Clusterization {
	const vector<FastqRead> &reads_;
	vector<pair<size_t, size_t> > index_cluster_;
    map<size_t, pair<string, size_t> > cluster_seq_size_;

public:
	Clusterization(const vector<FastqRead> &reads) :
		reads_(reads) { }

	void Add(size_t index, size_t cluster) {
        
		index_cluster_.push_back(make_pair(index, cluster));
	}

    void AddClusterSequence(size_t cluster, string seq, size_t size) {
        cluster_seq_size_[cluster].first = seq;
        cluster_seq_size_[cluster].second = size;
    }

	void WriteReadClusterMap(string fname) {
		ofstream out(fname.c_str());
		for(size_t i = 0; i < index_cluster_.size(); i++) {
            string name = reads_[index_cluster_[i].first].name.substr(1, reads_[index_cluster_[i].first].name.size() - 1);
			out << name << "\t" << index_cluster_[i].second << endl;
        }
	}

    void WriteClustersFasta(string fname) {
        ofstream out(fname.c_str());
        for(auto it = cluster_seq_size_.begin(); it != cluster_seq_size_.end(); it++) {
            out << ">cluster___" << it->first << "___size___" << it->second.second << endl;
            out << it->second.first << endl;
        }
    }

    size_t ClustersSize() { return cluster_seq_size_.size(); }
};

class Clusterizator {
	const vector<FastqRead> &reads_;
	ReadIndex<FastqRead> &read_index_;
	SequenceComparer &seq_comparer_;
	ConcurrentDSU clusters_;
    map<size_t, string> read_seq_map_;

    string GetSuperString(string s1, string s2, SeqComparisonResults comparison_result) {
        assert(comparison_result.match);
        if(comparison_result.second_reverse)
            s2 = reverse_complementary(s2);
        size_t overlap_size = comparison_result.s2_overlap_pos.second - comparison_result.s2_overlap_pos.first + 1;
        string superstring;

        // substring before overlap
        if(comparison_result.s1_overlap_pos.first != 0)
            superstring = s1.substr(0, comparison_result.s1_overlap_pos.first);
        else {
            if(comparison_result.s2_overlap_pos.first != 0)
                superstring = s2.substr(0, comparison_result.s2_overlap_pos.first);
        }

        // overlap substring
        superstring = superstring + s1.substr(comparison_result.s1_overlap_pos.first, overlap_size);

        // substring after overlap
        if(comparison_result.s1_overlap_pos.second != s1.size() - 1)
            superstring = superstring + s1.substr(comparison_result.s1_overlap_pos.second + 1, s1.size() - comparison_result.s1_overlap_pos.second - 1);
        else
            if(comparison_result.s2_overlap_pos.second != s2.size() - 1)
                superstring = superstring + s2.substr(comparison_result.s2_overlap_pos.second + 1, s2.size() - comparison_result.s2_overlap_pos.second - 1);

        return superstring;

        if(comparison_result.s1_overlap_pos.first >= comparison_result.s2_overlap_pos.first) {
            superstring = s1.substr(0, comparison_result.s1_overlap_pos.second + 1);
            superstring = superstring + s2.substr(comparison_result.s2_overlap_pos.second + 1, s2.size() - 
                comparison_result.s2_overlap_pos.second - 1);
        }
        else {
            superstring = s2.substr(0, comparison_result.s2_overlap_pos.second + 1);
            superstring = superstring + s1.substr(comparison_result.s1_overlap_pos.second + 1, s1.size() - 
                comparison_result.s1_overlap_pos.second - 1);        
        }
        return superstring;
    }

    void InitializeReadSequenceMap() {
        for(size_t i = 0; i < reads_.size(); i++) {
            read_seq_map_[i] = reads_[i].seq;
        }
    }

public:
	Clusterizator(const vector<FastqRead> &reads,
			ReadIndex<FastqRead> &read_index,
			SequenceComparer &seq_comparer) :
	reads_(reads),
	read_index_(read_index),
	seq_comparer_(seq_comparer),
	clusters_(reads.size()) {	}

	Clusterization CostructClusters() {
        InitializeReadSequenceMap();
#pragma omp parallel for
		for(size_t i = 0; i < reads_.size(); i++) {
			auto candidates = read_index_.GetCandidatesFor(i);
            for(auto it = candidates.begin(); it != candidates.end(); it++) {
                size_t cluster1 = clusters_.find_set(i);
                size_t cluster2 = clusters_.find_set(*it);
				if(cluster1 != cluster2) {
                    string read_seq1 = read_seq_map_[cluster1];
                    string read_seq2 = read_seq_map_[cluster2];
                    auto comparison_result = seq_comparer_.SequencesMatch(read_seq1, read_seq2);
					if(comparison_result.match) {
                        //cout << "Sequences match" << endl << "Read1: " << read_seq1 << endl << "Read2: " << read_seq2 << endl;
                        //cout << "Cluster1: " << clusters_.find_set(i) << " , cluster2: " << clusters_.find_set(*it) << endl;
						clusters_.unite(i, *it);
                        string superstring = GetSuperString(read_seq1, read_seq2, comparison_result);
                        //cout << "Superstring: " << superstring << endl;
                        size_t new_cluster = clusters_.find_set(i);
                        //cout << "New cluster: " << new_cluster << endl;
                        read_seq_map_[new_cluster] = superstring;
					}
				}
			}
		}
		cout << clusters_.num_sets() << " clusters were constructed" << endl;

		Clusterization result(reads_);
		for(size_t i = 0; i < reads_.size(); i++)
			result.Add(i, clusters_.find_set(i));

        for(auto it = read_seq_map_.begin(); it != read_seq_map_.end(); it++) {
            size_t cluster_id = clusters_.find_set(it->first);
            if(it->first == cluster_id) {
                result.AddClusterSequence(cluster_id, it->second, clusters_.set_size(cluster_id));
            }
        }

        assert(result.ClustersSize() == clusters_.num_sets());

		return result;

	}
};
