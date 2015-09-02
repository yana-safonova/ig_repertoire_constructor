from Bio import pairwise2

class GoodBarcodeMatch():
    def __init__(self, barcode_id, corr_cluster_id, all_corr_cluster_ids):
        self.barcode_id = barcode_id
        self.corr_cluster_id = corr_cluster_id
        self.all_corr_cluster_ids = all_corr_cluster_ids

class Distance():
    def __init__(self, mismatches, gaps):
        self.mismatches = mismatches
        self.gaps = gaps

class BarcodeMetrics():
    def __init__(self, barcode_repertoire, compressed_barcode_repertoire, barcode_cluster_matches, 
                       data_repertoire, compressed_data_repertoire, data_cluster_matches, rate_cutoff = 0.9):
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
        self.bad_barcodes = 0
        self.has_reads = bool(self.barcode_rep.read_clusters and self.data_rep.read_clusters)
        self.rate_cutoff = rate_cutoff
        self.good_barcode_matches= []

    def calculate_distances(self):
        distances = []
        '''
        handler = open("dist1.txt", 'w')
        handler.write("Barcode_id\tSize\tAbundance\tMatches\n")
        '''
        dist1_sizes = []
        dist1_abundances = []
        for cluster_id, cluster in self.compressed_barcode_rep.clusters.items():
            if cluster_id in self.barcode_cluster_matches:
                for i in range(cluster.abundance):
                    distances.append(self.barcode_cluster_matches[cluster_id][1])
                if self.barcode_cluster_matches[cluster_id][1] == 1:
                    dist1_sizes.append(cluster.size)
                    dist1_abundances.append(cluster.abundance)
        '''
                    handler.write(str(cluster_id) + '\t' + str(cluster.size) + '\t' + str(cluster.abundance) + '\t' +
                                  '\t'.join(str(i) for i in self.barcode_cluster_matches[cluster_id][0]) + '\n')
        handler.write("Min/avg/max sizes are: " + 
                       str(min(dist1_sizes)) + " / " +
                       str(float(sum(dist1_sizes)) / len(dist1_sizes)) + " / " +
                       str(max(dist1_sizes)) + "\n")
        handler.write("Min/avg/max abundances are: " + 
                       str(min(dist1_abundances)) + " / " +
                       str(float(sum(dist1_abundances)) / len(dist1_abundances)) + " / " +
                       str(max(dist1_abundances)) + "\n")
        handler.close()
        '''
        for dist in range(max(distances) + 1):
            cnt = distances.count(dist)
            if cnt:
                self.distances_dict[dist] = cnt

    def calculate_good_barcodes_number(self):
        if not self.has_reads:
            return
        for cluster_id, reads in self.barcode_rep.cluster_reads.items():
            data_corr_clusters = dict() # cluster_id -> number_of_shared_reads
            good_reads_num = 0
            for read_id in reads:
                if read_id not in self.data_rep.read_clusters:
                    continue
                good_reads_num += 1
                corr_data_cluster_id = self.data_rep.read_clusters[read_id]
                if corr_data_cluster_id not in data_corr_clusters:
                    data_corr_clusters[corr_data_cluster_id] = 0
                data_corr_clusters[corr_data_cluster_id] += 1
            if not data_corr_clusters:
                continue
            corr_data_cluster_id, corr_cluster_shared_reads_num = \
                max(data_corr_clusters.items(), key = lambda a: a[1])
            corr_cluster_size = len(self.data_rep.cluster_reads[corr_data_cluster_id])
            if corr_cluster_shared_reads_num >= self.rate_cutoff * good_reads_num:
                self.good_barcodes += 1
                match = GoodBarcodeMatch(cluster_id, corr_data_cluster_id, data_corr_clusters.keys())
                self.good_barcode_matches.append(match)
            else:
                self.bad_barcodes += 1

    def get_distance_for_alignment(self, alignment):
        mismatches = 0
        gaps = 0
        first_i = max(next(i for i, c in enumerate(alignment[0]) if c != '-'), 
                  next(i for i, c in enumerate(alignment[1]) if c != '-'))
        last_i = min(max(i for i, c in alignment[0] if c != '-'), 
                 max(i for i, c in alignment[1] if c != '-'))
        for i in len(alignment[0]):
            if i < first_i or i > last_i:
                continue
            if alignment[0][i] == '-' or alignment[1][i] == '-':
                gaps += 1
                continue
            if alignment[0][i] != alignment[1][i]:
                mismatches += 1
        return Distance(mismatches, gaps)

    def process_good_barcode_matches(self):
        distances = []
        for barcode_match in self.good_barcode_matches:
            for cluster_id in barcode_match.all_corr_cluster_ids:
                if cluster_id == barcode_match.corr_cluster_id:
                    continue
                alignment = pairwise2.align.globalxx(self.barcode_rep.clusters[barcode_match.barcode_id].seq,
                                         self.data_rep.clusters[cluster_id].seq, one_alignment_only=True)[0]
                distance = self.get_distance_for_alignment(alignment)
                distances.append(distance)
        self.min_gaps_distance = min(d.gaps for d in distances)
        self.avg_gaps_distance = float(sum(d.gaps for d in distances)) / len(distances)
        self.max_gaps_distance = max(d.gaps for d in distances)
        self.min_mismatches_distance = min(d.mismatches for d in distances)
        self.avg_mismatches_distance = float(sum(d.mismatches for d in distances)) / len(distances)
        self.max_mismatches_distance = max(d.mismatches for d in distances)

    def compute_lengths(self):
        barcode_lengths = self.barcode_rep.get_all_cluster_seq_lengths()
        self.min_barcode_length = min(barcode_lengths)
        self.avg_barcode_length = float(sum(barcode_lengths)) / len(barcode_lengths)
        self.max_barcode_length = max(barcode_lengths)
        data_clusters_lengths = self.data_rep.get_all_cluster_seq_lengths()
        self.min_data_clusters_length = min(data_clusters_lengths)
        self.avg_data_clusters_length = float(sum(data_clusters_lengths)) / len(data_clusters_lengths)
        self.max_data_clusters_length = max(data_clusters_lengths)

    def evaluate(self):
        self.calculate_distances()
        self.compute_lengths()
        if self.has_reads:
            self.calculate_good_barcodes_number()
            # self.process_good_barcode_matches()

    def write(self, filename):
        handler = open(filename, 'w')
        handler.write('Total number of barcodes: ' + str(self.barcodes_number) + '\n')
        handler.write('Min barcode sequence length: ' + str(self.min_barcode_length) + '\n')
        handler.write('Avg barcode sequence length: ' + str(self.avg_barcode_length) + '\n')
        handler.write('Max barcode sequence length: ' + str(self.max_barcode_length) + '\n')
        '''
        handler.write('Total number of barcodes in compressed repertoire: ' + \
                      str(self.compressed_barcodes_number) + '\n')
        '''
        handler.write('\n')
        handler.write('Total number of clusters in data: ' + str(self.data_clusters_number) + '\n')
        handler.write('Min data cluster sequence length: ' + str(self.min_data_clusters_length) + '\n')
        handler.write('Avg data cluster sequence length: ' + str(self.avg_data_clusters_length) + '\n')
        handler.write('Max data cluster sequence length: ' + str(self.max_data_clusters_length) + '\n')
        '''
        handler.write('Total number of clusters in compressed repertoire: ' + \
                      str(self.compressed_data_clusters_number) + '\n')
        '''
        handler.write('\n')
        for dist, count in self.distances_dict.items():
            handler.write('Number of barcodes reconstructed with distance ' + str(dist) + 
                          ': ' + str(count) + '(' + 
                          str(round(float(count) / self.barcodes_number * 100, 2)) + '%)\n')
        handler.write('Number of isolated clusters in barcode: ' + \
            str(len(self.compressed_barcode_rep.get_isolated_cluster_sizes(self.barcode_cluster_matches))) + '\n')
        handler.write('Min size of isolated cluster in barcode: ' + \
            str(self.compressed_barcode_rep.get_min_isolated_cluster_size(self.barcode_cluster_matches)) + '\n') 
        handler.write('Avg size of isolated cluster in barcode: ' + \
            str(self.compressed_barcode_rep.get_avg_isolated_cluster_size(self.barcode_cluster_matches)) + '\n') 
        handler.write('Max size of isolated cluster in barcode: ' + \
            str(self.compressed_barcode_rep.get_max_isolated_cluster_size(self.barcode_cluster_matches)) + '\n') 
        handler.write('Number of isolated clusters in data: ' + \
            str(len(self.compressed_data_rep.get_isolated_cluster_sizes(self.data_cluster_matches))) + '\n')
        handler.write('Min size of isolated cluster in data: ' + \
            str(self.compressed_data_rep.get_min_isolated_cluster_size(self.data_cluster_matches)) + '\n') 
        handler.write('Avg size of isolated cluster in data: ' + \
            str(self.compressed_data_rep.get_avg_isolated_cluster_size(self.data_cluster_matches)) + '\n') 
        handler.write('Max size of isolated cluster in data: ' + \
            str(self.compressed_data_rep.get_max_isolated_cluster_size(self.data_cluster_matches)) + '\n') 
        if self.has_reads:
            handler.write('Good barcodes number: ' + str(self.good_barcodes) + '\n')
            handler.write('Bad barcodes number: ' + str(self.bad_barcodes) + '\n')
            '''
            handler.write('Min gap distance to good barcodes: ' + str(self.min_gaps_distance) + '\n')
            handler.write('Avg gap distance to good barcodes: ' + str(self.min_gaps_distance) + '\n')
            handler.write('Max gap distance to good barcodes: ' + str(self.min_gaps_distance) + '\n')
            handler.write('Min mismatche distance to good barcodes: ' + str(self.min_mismatches_distance) + '\n')
            handler.write('Avg mismatche distance to good barcodes: ' + str(self.min_mismatches_distance) + '\n')
            handler.write('Max mismatche distance to good barcodes: ' + str(self.min_mismatches_distance) + '\n')
            '''
        handler.close()
