#pragma once

#include "clusters.hpp"

struct Metrics {

	struct OriginalClustersMetrics {
		size_t num_clusters;
		size_t num_not_merged;
		size_t num_nm_big_singletons;
		size_t num_singletons;
		size_t max_cluster;

		OriginalClustersMetrics() :
			num_clusters(),
			num_not_merged(),
			num_nm_big_singletons(),
			num_singletons(),
			max_cluster() { }

		void Print(ostream &out) {
			out << "#original clusters\t\t\t\t" << num_clusters << endl;
			out << "#not merged\t\t\t\t" << num_not_merged << endl;
			out << "#not merged (not trivial + singletons)\t\t\t\t" << num_nm_big_singletons << endl;
			out << "#singletons\t\t\t\t" << num_singletons << endl;
			out << "max cluster\t\t\t\t" << max_cluster << endl;
		}
	};

	struct ConstructedClustersMetrics {
		size_t num_clusters;
		size_t num_errors;
		size_t num_singletons;
		size_t num_corr_singletons;
		size_t max_cluster;

		ConstructedClustersMetrics() :
			num_clusters(),
			num_errors(),
			num_singletons(),
			num_corr_singletons(),
			max_cluster() { }

		void Print(ostream &out) {
			out << "#constructed clusters\t\t\t\t" << num_clusters << endl;
			out << "#errors\t\t\t\t" << num_errors << endl;
			out << "#singletons\t\t\t\t" << num_singletons << endl;
			out << "max cluster\t\t\t\t" << max_cluster << endl;
		}
	};

	ConstructedClustersMetrics constructed_metrics;
	OriginalClustersMetrics original_metrics;

	double avg_fillin;
	double max_cluster_fillin;
	size_t correct_singletons;
	double used_reads;
	size_t lost_clusters_number;
	double lost_clusters_size;

	Metrics() :
		original_metrics(),
		constructed_metrics(),
		avg_fillin(),
		max_cluster_fillin(),
		correct_singletons(),
		used_reads(),
		lost_clusters_number(),
		lost_clusters_size() { }

	void Print(ostream &out) {
		original_metrics.Print(out);
		constructed_metrics.Print(out);
		out << "avg fill-in\t\t\t\t" << avg_fillin << endl;
		out << "max cluster fill-in\t\t\t\t" << max_cluster_fillin<< endl;
		out << "#correct singletons\t\t\t\t" << correct_singletons << endl;
		out << "used reads (%)\t\t\t\t" << used_reads << endl;
		out << "#lost clusters\t\t\t\t" << lost_clusters_number << endl;
		out << "lost clusters size (%)\t\t\t\t" << lost_clusters_size << endl;
	}
};

class ClustersEvaluator {
	Clusters original_clusters_;
	Clusters constructed_clusters_;
	string constructed_fname_;

	vector<string> original_reads_;
	vector<string> constructed_reads_;

	Metrics metrics_;

	void ComputeOriginalSingletons() {
		for(auto it = original_clusters_.clusters_begin(); it != original_clusters_.clusters_end(); it++)
			if(original_clusters_.IsClusterSingleton(it->first))
				metrics_.original_metrics.num_singletons++;
	}

	void ComputeConstructedSingletons() {
		for(auto it = constructed_clusters_.clusters_begin(); it != constructed_clusters_.clusters_end(); it++)
					if(constructed_clusters_.IsClusterSingleton(it->first)) {
						metrics_.constructed_metrics.num_singletons++;
						auto read_name = constructed_reads_[(it->second)[0]];
						if(original_clusters_.IsClusterSingleton(original_clusters_.GetCluster(read_name)))
							metrics_.constructed_metrics.num_corr_singletons++;
					}
		metrics_.correct_singletons = metrics_.constructed_metrics.num_corr_singletons;
	}

	void ComputeSingletonMetrics() {
		ComputeOriginalSingletons();
		ComputeConstructedSingletons();
	}

	void ComputeMaxClusterMetrics() {
		size_t max_orig_cluster = 0;
		for(auto it = original_clusters_.clusters_begin(); it != original_clusters_.clusters_end(); it++)
			if(metrics_.original_metrics.max_cluster < original_clusters_.ClusterSize(it->first)) {
				metrics_.original_metrics.max_cluster = original_clusters_.ClusterSize(it->first);
				max_orig_cluster = it->first;
			}

		for(auto it = constructed_clusters_.clusters_begin(); it != constructed_clusters_.clusters_end(); it++) {
			metrics_.constructed_metrics.max_cluster = max<size_t>(metrics_.constructed_metrics.max_cluster,
					constructed_clusters_.ClusterSize(it->first));
		}

		metrics_.max_cluster_fillin = ComputeNotMergedFillin(max_orig_cluster);
	}

	void FillBaseMetrics() {
		metrics_.original_metrics.num_clusters = original_clusters_.ClustersNumber();
		metrics_.constructed_metrics.num_clusters = constructed_clusters_.ClustersNumber();
		ComputeSingletonMetrics();
		ComputeMaxClusterMetrics();
		metrics_.used_reads = double(constructed_reads_.size()) / original_reads_.size() * 100;
	}

	vector<string> GetReadNamesFromCluster(size_t cluster_ind, Clusters &clusters, vector<string> &reads) {
		vector<size_t> read_inds = clusters.GetReadsByCluster(cluster_ind);
		vector<string> read_names;
		for(auto it = read_inds.begin(); it != read_inds.end(); it++)
			read_names.push_back(reads[*it]);
		return read_names;
	}

	bool IsClusterNotMerged(size_t orig_cluster_ind) {
		vector<string> read_names = GetReadNamesFromCluster(orig_cluster_ind, original_clusters_, original_reads_);
		set<size_t> constr_cluster_inds;
		for(auto it = read_names.begin(); it != read_names.end(); it++)
			if(constructed_clusters_.ReadExists(*it))
				constr_cluster_inds.insert(constructed_clusters_.GetCluster(*it));
		return constr_cluster_inds.size() > 1;
	}

	bool IsClusterBigAndSingletons(size_t orig_cluster_ind) {
		assert(IsClusterNotMerged(orig_cluster_ind));
		vector<string> read_names = GetReadNamesFromCluster(orig_cluster_ind, original_clusters_, original_reads_);

		set<size_t> not_trivial_clusters;
		for(auto it = read_names.begin(); it != read_names.end(); it++) {
			if(constructed_clusters_.GetClusterSizeByReadName(*it) > 1) {
				not_trivial_clusters.insert(constructed_clusters_.GetCluster(*it));
			}
		}
		return not_trivial_clusters.size() == 1;
	}

	bool IsConstructedClusterNotErroneous(size_t constr_cluster_ind) {
		vector<string> read_names = GetReadNamesFromCluster(constr_cluster_ind,
				constructed_clusters_, constructed_reads_);
		set<size_t> orig_clusters_inds;
		for(auto it = read_names.begin(); it != read_names.end(); it++)
			orig_clusters_inds.insert(original_clusters_.GetCluster(*it));
		return orig_clusters_inds.size() == 1;
	}

	size_t GetOriginalSuperCluster(size_t constr_cluster_ind) {
		vector<string> read_names = GetReadNamesFromCluster(constr_cluster_ind,
				constructed_clusters_, constructed_reads_);
		set<size_t> orig_clusters_inds;
		for(auto it = read_names.begin(); it != read_names.end(); it++)
			orig_clusters_inds.insert(original_clusters_.GetCluster(*it));
		assert(orig_clusters_inds.size() == 1);
		return *(orig_clusters_inds.begin());
	}

	double ComputeNotMergedFillin(size_t orig_cluster_ind) {
		vector<string> read_names = GetReadNamesFromCluster(orig_cluster_ind, original_clusters_, original_reads_);
		map<size_t, vector<size_t> > constr_subclusters;
		for(auto it = read_names.begin(); it != read_names.end(); it++) {
			size_t cluster_ind = constructed_clusters_.GetCluster(*it);
			constr_subclusters[constructed_clusters_.ClusterSize(cluster_ind)].push_back(cluster_ind);
		}

		for(auto it = constr_subclusters.rbegin(); it != constr_subclusters.rend(); it++) {
			auto clusters = it->second;
			for(auto cluster_ind = clusters.begin(); cluster_ind != clusters.end(); cluster_ind++)
				if(IsConstructedClusterNotErroneous(*cluster_ind)) {
					return double(constructed_clusters_.ClusterSize(*cluster_ind)) / original_clusters_.ClusterSize(orig_cluster_ind);
				}
		}
		return 0;
	}

	void ComputeNotMergedMetrics() {
		double sum_fill = 0;

		for(auto it = original_clusters_.clusters_begin(); it != original_clusters_.clusters_end(); it++) {
			size_t orig_cluster_ind = it->first;
			if(original_clusters_.IsClusterSingleton(orig_cluster_ind))
				continue;

			// for each non trivial cluster we compute it "not merged"
			if(!IsClusterNotMerged(orig_cluster_ind))
				continue;

			metrics_.original_metrics.num_not_merged++;
			if(IsClusterBigAndSingletons(orig_cluster_ind))
				metrics_.original_metrics.num_nm_big_singletons++;

			size_t orig_cluster_size = original_clusters_.ClusterSize(orig_cluster_ind);
			double orig_cluster_fill = ComputeNotMergedFillin(orig_cluster_ind);

			sum_fill += orig_cluster_fill;
		}
		metrics_.avg_fillin = sum_fill / metrics_.original_metrics.num_not_merged;
	}

	void ComputeErrorMetrics() {
		metrics_.constructed_metrics.num_errors = 0;
		for(auto it = constructed_clusters_.clusters_begin(); it != constructed_clusters_.clusters_end(); it++) {
			if(!IsConstructedClusterNotErroneous(it->first))
				metrics_.constructed_metrics.num_errors++;
			else {
				size_t orig_supercluster = GetOriginalSuperCluster(it->first);
				if(constructed_clusters_.ClusterSize(it->first) > original_clusters_.ClusterSize(orig_supercluster))
				metrics_.constructed_metrics.num_errors++;
			}
		}
	}

	bool IsOriginalClusterLost(size_t orig_cluster_ind) {
		auto reads = original_clusters_.GetReadsByCluster(orig_cluster_ind);
		for(auto read = reads.begin(); read != reads.end(); read++)
			if(constructed_clusters_.ReadExists(original_reads_[*read]))
				return false;
		return true;
	}

	void ComputeLostClustersMetrics() {
		size_t lost_cluster_size = 0;
		for(auto it = original_clusters_.clusters_begin(); it != original_clusters_.clusters_end();
				it++) {
			if(IsOriginalClusterLost(it->first)) {
				metrics_.lost_clusters_number++;
				lost_cluster_size += original_clusters_.ClusterSize(it->first);
			}
		}
		metrics_.lost_clusters_size = double(lost_cluster_size) / original_reads_.size();
	}

public:
	ClustersEvaluator(char *original_fname, char *constructed_fname) :
		original_clusters_(original_fname),
		constructed_clusters_(constructed_fname),
		constructed_fname_(constructed_fname) {
		cout << "Extraction of original clusters from " << original_fname << endl;
		original_clusters_.ExtractFromFile();
		cout << "Extraction of constructed clusters from " << constructed_fname << endl;
		constructed_clusters_.ExtractFromFile();
	}

	void Evaluate() {
		original_reads_ = original_clusters_.Reads();
		constructed_reads_ = constructed_clusters_.Reads();

		FillBaseMetrics();
		ComputeNotMergedMetrics();
		ComputeErrorMetrics();
		ComputeLostClustersMetrics();
	}

	Metrics GetMetrics() { return metrics_; }

};
