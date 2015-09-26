#pragma once

#include "../include_me.hpp"
#include "../string_tools.hpp"

class Clusters {
	ifstream src_;

	map<string, size_t> read_cluster_;
	map<string, size_t> read_index_;
	vector<string> reads_;
	map<size_t, vector<size_t> > clusters_;

public:
	Clusters(string fname) :
		src_(fname.c_str()) {
		assert(!src_.fail());
	}

	void ExtractFromFile() {
		while(!src_.eof()) {
			string tmp;
			getline(src_, tmp);
			if(tmp != "") {
				vector<string> splits = split(tmp, '\t');
                if(splits.size() != 2) {
                    cout << "String " << tmp << " is wrong. Nummber of splits is " << splits.size() << endl;
				    assert(splits.size() == 2);
                }
				string name;
				size_t pos = splits[0].find("fr=");
				name = splits[0].substr(0, pos);
				read_cluster_[/*splits[0]*/name] = string_to_number<size_t>(splits[1]);
			}
		}

		size_t index = 0;
		for(auto it = read_cluster_.begin(); it != read_cluster_.end(); it++) {
			read_index_[it->first] = index;
			index++;
			reads_.push_back(it->first);
		}

		for(auto it = read_cluster_.begin(); it != read_cluster_.end(); it++) {
			clusters_[it->second].push_back(read_index_[it->first]);
		}
	}

	vector<string> Reads() { return reads_; }

	map<size_t, vector<size_t> >::iterator clusters_begin() { return clusters_.begin(); }

	map<size_t, vector<size_t> >::iterator clusters_end() { return clusters_.end(); }

	size_t GetCluster(size_t read_index) {
		assert(read_index < reads_.size());
		return read_cluster_[reads_[read_index]];
	}

	size_t GetCluster(string read_name) {
		if(read_index_.find(read_name) == read_index_.end())
			return size_t(-1);
		return GetCluster(read_index_[read_name]);
	}

	size_t GetClusterSizeByReadName(string read_name) {
		if(read_index_.find(read_name) == read_index_.end())
			return size_t(-1);
		return ClusterSize(GetCluster(read_name));
	}

	size_t ClusterSize(size_t cluster) {
		return clusters_[cluster].size();
	}

	vector<size_t> GetReadsByCluster(size_t cluster) {
		return clusters_[cluster];
	}

	size_t ClustersNumber() {
		return clusters_.size();
	}

	bool IsClusterSingleton(size_t index) {
		return clusters_[index].size() == 1;
	}

	bool ReadExists(string read_name) {
		return read_index_.find(read_name) != read_index_.end();
	}

};
