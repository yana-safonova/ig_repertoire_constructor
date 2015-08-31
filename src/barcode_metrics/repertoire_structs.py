import os.path 

class Repertoire:
    def __init__(self, clusters_filename, rcm_filename, clusters, cluster_reads = None, read_clusters = None, name = None):
        self.clusters_filename = clusters_filename
        self.rcm_filename = rcm_filename
        self.clusters = clusters
        self.cluster_reads = cluster_reads
        self.read_clusters = read_clusters
        self.name = name if name else os.path.basename(clusters_filename)

    def __str__(self):
        return ('Clusters: ' + str(self.clusters) + ', \n reads: ' + 
                ('None' if not self.cluster_reads else self.cluster_reads))

    def __repr__(self):
        return str(self)

    def get_cluster_seq_length(self, cluster_id):
        return len(self.clusters[cluster_id].seq)

    def get_isolated_cluster_sizes(self, cluster_matches):
        sizes = []
        for cluster_id, cluster in self.clusters.items():
            if cluster_id not in cluster_matches:
                sizes.append(cluster.size)
        return sizes

    def get_min_isolated_cluster_size(self, cluster_matches):
        sizes = self.get_isolated_cluster_sizes(cluster_matches)
        if not sizes:
            return 0
        return min(self.get_isolated_cluster_sizes(cluster_matches))

    def get_max_isolated_cluster_size(self, cluster_matches):
        sizes = self.get_isolated_cluster_sizes(cluster_matches)
        if not sizes:
            return 0
        return max(self.get_isolated_cluster_sizes(cluster_matches))

    def get_avg_isolated_cluster_size(self, cluster_matches):
        sizes = self.get_isolated_cluster_sizes(cluster_matches)
        if not sizes:
            return 0
        return float(sum(sizes)) / len(sizes)

class Cluster:
    def __init__(self, size, seq, abundance = 1):
        self.size = size
        self.seq = seq
        self.abundance = abundance

    def __str__(self):
        return 'size ' + str(self.size) + ', seq ' + str(self.seq)

    def __repr__(self):
        return str(self)
