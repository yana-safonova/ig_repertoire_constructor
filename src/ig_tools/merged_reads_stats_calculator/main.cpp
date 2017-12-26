#include "omp.h"
#include "quality_statistics.hpp"

void OutputMergedRL(std::vector<FastqRead> &reads, std::ostream &out) {
	for(auto it = reads.begin(); it != reads.end(); it++)
		out << it->seq.size() << std::endl;
}

int main(int argc, char *argv[]) {
	omp_set_num_threads(1);

	/*
	 * argv[1] - left raw reads
	 * argv[2] - right raw reads
	 * argv[3] - merged reads
	 */

	if(argc != 5) {
		std::cout << "Invalid input parameters" << std::endl <<
				"\t./compute_quality_stats left_reads.fq right_reads.fq merged_reads.fq output_dir" << std::endl;
		return 1;
	}

	size_t phred_offset = 33;

	//cout << "Statistics for paired reads" << endl;
	auto paired_reads = PairedFastqReader(argv[1], argv[2]).Read();
	PairedReadQialityStatsCalculator paired_calculator(paired_reads, phred_offset);
	paired_calculator.Calculate();
	//paired_calculator.Stats().ShortPrint(cout);
	std::ofstream out1((std::string(argv[4]) + "/paired_nucl_qual.stats").c_str());
	paired_calculator.Stats().OutputNucleotideQuality(out1);

	//cout << "Statistics for paired reads" << endl;
	auto merged_reads = SingleFastqReader(argv[3]).ReadFile();
	MergedReadQualityStatsCalculator merged_calculator(merged_reads, phred_offset);
	merged_calculator.Calculate();
	//merged_calculator.Stats().ShortPrint(cout);
	std::ofstream out2((std::string(argv[4]) + "/merged_nucl_qual.stats").c_str());
	merged_calculator.Stats().OutputNucleotideQuality(out2);

	std::ofstream out3((std::string(argv[4]) + "/merged_rl.stats").c_str());
	OutputMergedRL(merged_reads, out3);

	return 0;
}
