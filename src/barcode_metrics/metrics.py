from Bio import pairwise2
from scipy.stats.stats import pearsonr

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
                       data_repertoire, compressed_data_repertoire, data_cluster_matches, rate_cutoff):
        self.barcode_rep = barcode_repertoire
        self.compressed_barcode_rep = compressed_barcode_repertoire
        self.barcode_cluster_matches = barcode_cluster_matches
        self.data_rep = data_repertoire
        self.compressed_data_rep = compressed_data_repertoire
        self.data_cluster_matches = data_cluster_matches
        self.rate_cutoff = rate_cutoff
        self.distances_dict = {} # distance -> count
        self.barcodes_number = len(self.barcode_rep.clusters)
        self.data_clusters_number = len(self.data_rep.clusters)
        self.has_reads = bool(self.barcode_rep.read_clusters and self.data_rep.read_clusters)
        self.good_barcode_matches= []
        self.has_barcodes_dist = False
        # self.corresponding_cluster_sizes = ([], [])
        # self.sizes_correlation = 0
        self.sep_metric_names = ["#clusters", "#isolated clusters",
            "#clusters (>=10)", "#clusters (>=50)", "#clusters (>=100)", "#clusters (>=500)", "#clusters (>=1000)",
            "Min sequence length", "Avg sequence length", "Max sequence length",
            "Min size of cluster", "Avg size of cluster", "Max size of cluster",
            "Min size of isolated cluster", "Avg size of isolated cluster", "Max size of isolated cluster"]
        self.sep_metrics_dict = {}
        self.gen_metrics_dict = {}

    def calculate_distances(self, log):
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
                if self.barcode_cluster_matches[cluster_id][1] == 0:
                    corr_size = sum(self.compressed_data_rep.clusters[c].size for c in
                                    self.barcode_cluster_matches[cluster_id][0])
                '''
                    self.corresponding_cluster_sizes[0].append(self.compressed_barcode_rep.clusters[cluster_id].size)
                    self.corresponding_cluster_sizes[1].append(corr_size)
                if self.barcode_cluster_matches[cluster_id][1] == 1:
                    dist1_sizes.append(cluster.size)
                    dist1_abundances.append(cluster.abundance)
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
        # self.sizes_correlation = pearsonr(self.corresponding_cluster_sizes[0], self.corresponding_cluster_sizes[1])[0]
        for dist in range(max(distances) + 1):
            cnt = distances.count(dist)
            if cnt:
                self.distances_dict[dist] = cnt

    def calculate_good_barcodes_number(self, log):
        good_barcodes = 0
        bad_barcodes = 0
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
                good_barcodes += 1
                match = GoodBarcodeMatch(cluster_id, corr_data_cluster_id, data_corr_clusters.keys())
                self.good_barcode_matches.append(match)
            else:
                bad_barcodes += 1
        self.gen_metrics_dict["#good barcodes"] = "{:d} ({:.2f}%)".format(good_barcodes, 100.0 * good_barcodes / self.barcodes_number)
        self.gen_metrics_dict["#bad barcodes"] = "{:d} ({:.2f}%)".format(bad_barcodes, 100.0 * bad_barcodes/ self.barcodes_number)

    def get_distance_for_alignment(self, barcode_seq, data_seq):
        mismatches = 0
        gaps = 0
        first_i = max(next(i for i, c in enumerate(barcode_seq) if c != '-'), 
                  next(i for i, c in enumerate(data_seq) if c != '-'))
        last_i = min(max(i for i, c in enumerate(barcode_seq) if c != '-'), 
                 max(i for i, c in enumerate(data_seq) if c != '-'))
        for i in range(len(barcode_seq)):
            if i < first_i or i > last_i:
                continue
            if barcode_seq[i] == '-' or data_seq[i] == '-':
                gaps += 1
                continue
            if barcode_seq[i] != data_seq[i]:
                mismatches += 1
        return Distance(mismatches, gaps)

    def get_distances_for_alignment(self, alignment):
        distances = []
        barcode_record = alignment._records[0]
        for record in alignment._records[1:]:
            distances.append(self.get_distance_for_alignment(barcode_record.seq, record.seq))
        return distances

    '''
    def align_with_muscle(self, string_to_align):
        from Bio.Align.Applications import MuscleCommandline
        from Bio import AlignIO
        from StringIO import StringIO
        muscle_cline = MuscleCommandline(clwstrict=True, diags=True, maxiters=2)
        handler = StringIO()
        handler.write(string_to_align)
        stdout, stderr = muscle_cline(stdin=handler.getvalue())
        return AlignIO.read(StringIO(stdout), 'clustal')
    '''

    def process_good_barcode_matches(self, log):
        distances = []
        log.info("..Number of good barcodes is " + str(len(self.good_barcode_matches)))
        for i, barcode_match in enumerate(self.good_barcode_matches):
            if i % 100 == 0:
                log.info("..Processed " + str(i) + " good barcodes")
            for cluster_id in barcode_match.all_corr_cluster_ids:
                if cluster_id == barcode_match.corr_cluster_id:
                    continue
                if not self.barcode_rep.clusters[barcode_match.barcode_id].seq:
                    continue
                alignment = pairwise2.align.globalxx(self.barcode_rep.clusters[barcode_match.barcode_id].seq,
                    self.data_rep.clusters[cluster_id].seq, one_alignment_only=True)[0]
                distance = self.get_distance_for_alignment(alignment[0], alignment[1])
                distances.append(distance)
        self.gen_metrics_dict["Min gap distance to good barcodes"] = \
            min(d.gaps for d in distances)
        self.gen_metrics_dict["Avg gap distance to good barcodes"] = \
            float(sum(d.gaps for d in distances)) / len(distances)
        self.gen_metrics_dict["Max gap distance to good barcodes"] = \
            max(d.gaps for d in distances)
        self.gen_metrics_dict["Min mismatches distance to good barcodes"] = \
            min(d.mismatches for d in distances)
        self.gen_metrics_dict["Avg mismatches distance to good barcodes"] = \
            float(sum(d.mismatches for d in distances)) / len(distances)
        self.gen_metrics_dict["Max mismatches distance to good barcodes"] = \
            max(d.mismatches for d in distances)

    def compute_length_metrics(self):
        self.barcode_lengths = self.barcode_rep.get_all_cluster_seq_lengths()
        self.sep_metrics_dict["barcodes"]["Min sequence length"] = min(self.barcode_lengths)
        self.sep_metrics_dict["barcodes"]["Avg sequence length"] = \
            str(round(float(sum(self.barcode_lengths)) / len(self.barcode_lengths), 2))
        self.sep_metrics_dict["barcodes"]["Max sequence length"] = max(self.barcode_lengths)
        self.data_clusters_lengths = self.data_rep.get_all_cluster_seq_lengths()
        self.sep_metrics_dict["data"]["Min sequence length"] = min(self.data_clusters_lengths)
        self.sep_metrics_dict["data"]["Avg sequence length"] = \
            str(round(float(sum(self.data_clusters_lengths)) / len(self.data_clusters_lengths), 2))
        self.sep_metrics_dict["data"]["Max sequence length"] = max(self.data_clusters_lengths)

    def compute_clusters_number_metrics(self):
        self.sep_metrics_dict["barcodes"]["#clusters"] = self.barcodes_number
        self.sep_metrics_dict["data"]["#clusters"] = self.data_clusters_number
        self.sep_metrics_dict["barcodes"]["#clusters (>=10)"] = len([s for s in 
            self.barcode_rep.get_all_cluster_sizes() if s >= 10])
        self.sep_metrics_dict["barcodes"]["#clusters (>=50)"] = len([s for s in 
            self.barcode_rep.get_all_cluster_sizes() if s >= 50])
        self.sep_metrics_dict["barcodes"]["#clusters (>=100)"] = len([s for s in 
            self.barcode_rep.get_all_cluster_sizes() if s >= 100])
        self.sep_metrics_dict["barcodes"]["#clusters (>=500)"] = len([s for s in 
            self.barcode_rep.get_all_cluster_sizes() if s >= 500])
        self.sep_metrics_dict["barcodes"]["#clusters (>=1000)"] = len([s for s in 
            self.barcode_rep.get_all_cluster_sizes() if s >= 1000])
        self.sep_metrics_dict["data"]["#clusters (>=10)"] = len([s for s in 
            self.data_rep.get_all_cluster_sizes() if s >= 10])
        self.sep_metrics_dict["data"]["#clusters (>=50)"] = len([s for s in 
            self.data_rep.get_all_cluster_sizes() if s >= 50])
        self.sep_metrics_dict["data"]["#clusters (>=100)"] = len([s for s in 
            self.data_rep.get_all_cluster_sizes() if s >= 100])
        self.sep_metrics_dict["data"]["#clusters (>=500)"] = len([s for s in 
            self.data_rep.get_all_cluster_sizes() if s >= 500])
        self.sep_metrics_dict["data"]["#clusters (>=1000)"] = len([s for s in 
            self.data_rep.get_all_cluster_sizes() if s >= 1000])
        isolated_barcodes = len(self.compressed_barcode_rep.get_isolated_cluster_sizes(self.barcode_cluster_matches))
        self.sep_metrics_dict["barcodes"]["#isolated clusters"] = "{:d} ({:.2f}%)".format(
            isolated_barcodes, isolated_barcodes / self.barcodes_number * 100.0)
        isolated_data_clusters = len(self.compressed_data_rep.get_isolated_cluster_sizes(self.data_cluster_matches))
        self.sep_metrics_dict["data"]["#isolated clusters"] = "{:d} ({:.2f}%)".format(
            isolated_data_clusters, isolated_data_clusters / self.data_clusters_number * 100.0)

    def compute_cluster_size_metrics(self):
        self.sep_metrics_dict["barcodes"]["Min size of cluster"] = \
            self.compressed_barcode_rep.get_min_cluster_size()
        self.sep_metrics_dict["data"]["Min size of cluster"] = \
            self.compressed_data_rep.get_min_cluster_size()
        self.sep_metrics_dict["barcodes"]["Avg size of cluster"] = \
            round(self.compressed_barcode_rep.get_avg_cluster_size(), 3)
        self.sep_metrics_dict["data"]["Avg size of cluster"] = \
            round(self.compressed_data_rep.get_avg_cluster_size(), 3)
        self.sep_metrics_dict["barcodes"]["Max size of cluster"] = \
            self.compressed_barcode_rep.get_max_cluster_size()
        self.sep_metrics_dict["data"]["Max size of cluster"] = \
            self.compressed_data_rep.get_max_cluster_size()
        self.sep_metrics_dict["barcodes"]["Min size of isolated cluster"] = \
            self.compressed_barcode_rep.get_min_isolated_cluster_size(self.barcode_cluster_matches)
        self.sep_metrics_dict["data"]["Min size of isolated cluster"] = \
            self.compressed_data_rep.get_min_isolated_cluster_size(self.data_cluster_matches)
        self.sep_metrics_dict["barcodes"]["Avg size of isolated cluster"] = \
            round(self.compressed_barcode_rep.get_avg_isolated_cluster_size(self.barcode_cluster_matches), 3)
        self.sep_metrics_dict["data"]["Avg size of isolated cluster"] = \
            round(self.compressed_data_rep.get_avg_isolated_cluster_size(self.data_cluster_matches), 3)
        self.sep_metrics_dict["barcodes"]["Max size of isolated cluster"] = \
            self.compressed_barcode_rep.get_max_isolated_cluster_size(self.barcode_cluster_matches)
        self.sep_metrics_dict["data"]["Max size of isolated cluster"] = \
            self.compressed_data_rep.get_max_isolated_cluster_size(self.data_cluster_matches)

    def compute_metrics(self):
        self.compute_clusters_number_metrics()
        self.compute_length_metrics()
        self.compute_cluster_size_metrics()

    def evaluate(self, log, compute_rcm_matches):
        self.sep_metrics_dict["barcodes"] = {}
        self.sep_metrics_dict["data"] = {}
        self.calculate_distances(log)
        self.compute_metrics()
        if self.has_reads:
            self.calculate_good_barcodes_number(log)
            if compute_rcm_matches:
                self.process_good_barcode_matches(log)

    def write(self, filename):
        handler = open(filename, 'w')
        handler.write("Barcodes -> " + self.barcode_rep.clusters_filename + "\n")
        handler.write("Data -> " + self.data_rep.clusters_filename + "\n")
        handler.write("\n")

        handler.write("{:<45s}{:<14s}{:s}\n".format("Metric name", "Barcodes", "Data"))
        fmt = "{:<45s}{:14s}{:s}\n"
        for metric_name in self.sep_metric_names:
            if metric_name not in self.sep_metrics_dict["barcodes"]:
                continue
            handler.write(fmt.format(metric_name, 
                str(self.sep_metrics_dict["barcodes"][metric_name]), 
                str(self.sep_metrics_dict["data"][metric_name])))

        if self.gen_metrics_dict:
            handler.write("\n{:<45s}{:s}\n".format("Metric name", "Value"))
        fmt = "{:<45s}{:s}\n"
        for metric_name in self.gen_metrics_dict:
            handler.write(fmt.format(metric_name, str(self.gen_metrics_dict[metric_name])))
        handler.close()

    '''
    def draw_sizes_correlation_plot(self, filename):
        import matplotlib.pyplot as plt
        plt.plot(self.corresponding_cluster_sizes[0], self.corresponding_cluster_sizes[1], 'b*')
        plt.savefig(filename)
    '''
