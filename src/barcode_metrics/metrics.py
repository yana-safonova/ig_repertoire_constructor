class BarcodeMetrics():
    def __init__(self, barcode_repertoire, compressed_barcode_repertoire, barcode_cluster_matches, 
                       data_repertoire, compressed_data_repertoire, data_cluster_matches):
        self.barcode_rep = barcode_repertoire
        self.compressed_barcode_rep = compressed_barcode_repertoire
        self.barcode_cluster_matches = barcode_cluster_matches
        self.data_rep = data_repertoire
        self.compressed_data_rep = compressed_data_repertoire
        self.data_cluster_matches = data_cluster_matches
        self.distances_dict = {} # distance -> count
        self.barcodes_number = len(self.barcode_rep.clusters)
        self.compressed_barcodes_number = len(self.compressed_barcode_rep.clusters)
        self.data_clusters_number = len(self.data_rep.clusters)
        self.compressed_data_clusters_number = len(self.compressed_data_rep.clusters)
        self.good_barcodes = 0
        self.erroneous_clusters = 0
        self.has_reads = bool(self.barcode_rep.read_clusters and self.data_rep.read_clusters)

    def calculate_distances(self):
        distances = []
        for cluster_id, cluster in self.compressed_barcode_rep.clusters.items():
            if cluster_id in self.barcode_cluster_matches:
                for i in range(cluster.abundance):
                    distances.append(self.barcode_cluster_matches[cluster_id][1])
        for dist in range(max(distances) + 1):
            cnt = distances.count(dist)
            if cnt:
                self.distances_dict[dist] = cnt

    def calculate_good_barcodes_number(self):
        if not self.has_reads:
            return
        for cluster_id, reads in self.barcode_rep.cluster_reads.items():
            data_corr_clusters = dict() # cluster_id -> number_of_shared_reads
            for read_id in reads:
                if read_id not in self.data_rep.read_clusters:
                    continue
                corr_data_cluster_id = self.data_rep.read_clusters[read_id]
                if corr_data_cluster_id not in data_corr_clusters:
                    data_corr_clusters[corr_data_cluster_id] = 0
                data_corr_clusters[corr_data_cluster_id] += 1
            if not data_corr_clusters:
                continue
            corr_data_cluster_id, corr_cluster_shared_reads_num = \
                max(data_corr_clusters.items(), key = lambda a: a[1])
            if corr_cluster_shared_reads_num >= 0.9 * len(reads):
                self.good_barcodes += 1

    def evaluate(self):
        self.calculate_distances()
        self.calculate_good_barcodes_number()

    def write(self, filename):
        handler = open(filename, 'w')
        handler.write('Total number of barcodes is ' + str(self.barcodes_number) + '\n')
        handler.write('Total number of barcodes in compressed repertoire is ' + \
                      str(self.compressed_barcodes_number) + '\n')
        handler.write('Total number of clusters in data is ' + str(self.data_clusters_number) + '\n')
        handler.write('Total number of clusters in compressed repertoire is ' + \
                      str(self.compressed_data_clusters_number) + '\n')
        handler.write('\n')
        for dist, count in self.distances_dict.items():
            handler.write('Number of barcodes reconstructed with distance ' + str(dist) + 
                          ' is ' + str(count) + '(' + 
                          str(round(float(count) / self.barcodes_number * 100, 2)) + '%)\n')
        handler.write('Number of isolated clusters in barcode is ' + \
            str(len(self.compressed_barcode_rep.get_isolated_cluster_sizes(self.barcode_cluster_matches))) + '\n')
        handler.write('Min size of isolated cluster in barcode is ' + \
            str(self.compressed_barcode_rep.get_min_isolated_cluster_size(self.barcode_cluster_matches)) + '\n') 
        handler.write('Avg size of isolated cluster in barcode is ' + \
            str(self.compressed_barcode_rep.get_avg_isolated_cluster_size(self.barcode_cluster_matches)) + '\n') 
        handler.write('Max size of isolated cluster in barcode is ' + \
            str(self.compressed_barcode_rep.get_max_isolated_cluster_size(self.barcode_cluster_matches)) + '\n') 
        handler.write('Number of isolated clusters in data is ' + \
            str(len(self.compressed_data_rep.get_isolated_cluster_sizes(self.data_cluster_matches))) + '\n')
        handler.write('Min size of isolated cluster in data is ' + \
            str(self.compressed_data_rep.get_min_isolated_cluster_size(self.data_cluster_matches)) + '\n') 
        handler.write('Avg size of isolated cluster in data is ' + \
            str(self.compressed_data_rep.get_avg_isolated_cluster_size(self.data_cluster_matches)) + '\n') 
        handler.write('Max size of isolated cluster in data is ' + \
            str(self.compressed_data_rep.get_max_isolated_cluster_size(self.data_cluster_matches)) + '\n') 
        if self.has_reads:
            handler.write('Good barcodes number is ' + str(self.good_barcodes) + '\n')
        handler.close()
