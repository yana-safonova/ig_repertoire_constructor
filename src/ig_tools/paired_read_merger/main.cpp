#include "reads_merger.hpp"

merger_setting parse_settings(int argc, char *argv[]) {
	merger_setting setting;
	string min_overlap_str = "--min-overlap=";
	string max_mismatch_str = "--max-mismatch=";
	string simulated_mode_str = "--simulated-mode";
	for(size_t i = 4; i < argc; i++) {
		string tmp(argv[i]);
		if(tmp.substr(0, min_overlap_str.size()) == min_overlap_str) {
			tmp = tmp.substr(min_overlap_str.size(), tmp.size() - min_overlap_str.size());
			setting.min_overlap = string_to_number<size_t>(tmp);
		}
		else if(tmp.substr(0, max_mismatch_str.size()) == max_mismatch_str) {
			tmp = tmp.substr(max_mismatch_str.size(), tmp.size() - max_mismatch_str.size());
			setting.max_mismatch_rate = string_to_number<double>(tmp);
		}
		else if(tmp == simulated_mode_str)
			setting.simulated_mode = true;

	}
	setting.print(cout);
	return setting;
}

int main(int argc, char *argv[]) {
	/*
	 * argv[1] - left fastq reads
	 * argv[2] - right fastq reads
	 * argv[3] - prefix of output files (prefix.fastq prefix.stats)
	 */
	if(argc < 4) {
		cout << "paired_read_merger left_reads.fq right_reads.fq output_prefix [--min-overlap=N1 --max-mismatch=N2 --simulated-mode]" << endl;
		return 1;
	}

	vector<paired_fastq_read> paired_reads = PairedFastqReader(argv[1],
			argv[2]).Read();
	cout << paired_reads.size() << " paired reads were read from " << argv[1] <<
			" and " << argv[2] << endl;
	vector<fastq_read> merged_reads = PairedReadsMerger(parse_settings(argc, argv)).Merge(paired_reads);
	cout << merged_reads.size() << " read from " << paired_reads.size() << " were successfully merged" << endl;
	FastqWriter(string(argv[3]) + ".fastq").Write(merged_reads);
	cout << "Merged reads were written to " << string(argv[3]) + ".fastq" << endl; 
	return 0;
}
