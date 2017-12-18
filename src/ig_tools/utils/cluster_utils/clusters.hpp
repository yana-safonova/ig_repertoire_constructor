#pragma once

#include <fstream>
#include <set>
#include <string>
#include "../string_tools.hpp"

class Clusters {
	std::ifstream src_;

	std::map<std::string, size_t> read_cluster_;
	std::map<std::string, size_t> read_index_;
	std::vector<std::string> reads_;
	std::map<size_t, std::vector<size_t> > clusters_;

public:
	Clusters(std::string fname) :
		src_(fname.c_str()) {
		VERIFY(!src_.fail());
	}

	void ExtractFromFile() {
		while(!src_.eof()) {
			std::string tmp;
			getline(src_, tmp);
			if(tmp != "") {
				std::vector<std::string> splits = split(tmp, '\t');
                if(splits.size() != 2) {
                    cout << "String " << tmp << " is wrong. Nummber of splits is " << splits.size() << endl;
					VERIFY(splits.size() == 2);
                }
				std::string name;
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

	std::vector<std::string> Reads() { return reads_; }

	std::map<size_t, std::vector<size_t> >::iterator clusters_begin() { return clusters_.begin(); }

	std::map<size_t, std::vector<size_t> >::iterator clusters_end() { return clusters_.end(); }

	size_t GetCluster(size_t read_index) {
		VERIFY(read_index < reads_.size());
		return read_cluster_[reads_[read_index]];
	}

	size_t GetCluster(std::string read_name) {
		if(read_index_.find(read_name) == read_index_.end())
			return size_t(-1);
		return GetCluster(read_index_[read_name]);
	}

	size_t GetClusterSizeByReadName(std::string read_name) {
		if(read_index_.find(read_name) == read_index_.end())
			return size_t(-1);
		return ClusterSize(GetCluster(read_name));
	}

	size_t ClusterSize(size_t cluster) {
		return clusters_[cluster].size();
	}

	std::vector<size_t> GetReadsByCluster(size_t cluster) {
		return clusters_[cluster];
	}

	size_t ClustersNumber() {
		return clusters_.size();
	}

	bool IsClusterSingleton(size_t index) {
		return clusters_[index].size() == 1;
	}

	bool ReadExists(std::string read_name) {
		return read_index_.find(read_name) != read_index_.end();
	}

};
