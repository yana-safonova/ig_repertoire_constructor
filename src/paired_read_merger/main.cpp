#include "reads_merger.hpp"

#include <boost/program_options.hpp>

namespace po = boost::program_options;
namespace rm = reads_merger;

void create_console_logger(string log_filename) {
	using namespace logging;
	logger *lg = create_logger(path::FileExists(log_filename) ? log_filename : "");
	lg->add_writer(std::make_shared<console_writer>());
	attach_logger(lg);
}

int main(int argc, char *argv[]) {

	vector<string> in_left_files, in_right_files;
	string out_good_file;
	string out_bad_file;
	rm::merger_setting settings;

	create_console_logger("");

	//parsing command line
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "Help message")
		("min-overlap", po::value<size_t>(&settings.min_overlap)->default_value(50), "set minimal overlap size")
		("max-mismatch", po::value<double>(&settings.max_mismatch_rate)->default_value(0.1), "set maximal mismatch rate")
		("1", po::value<vector<string>>(&in_left_files), "input files containing left reads")
		("2", po::value<vector<string>>(&in_right_files), "input files containing right reads")
		("o", po::value<string>(&out_good_file), "output file for merged reads")
		("bad", po::value<string>(&out_bad_file), "output file for non-merged reads")
		("qual", po::value<size_t>(&settings.base_quality)->default_value(33), "set base quality")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);    

	if (vm.count("help")) {
		cout << desc << "\n";
		return 1;
	}
	if (in_left_files.size() != in_right_files.size()) {
		INFO("Numbers of files with left and right reads are different.");
		return 2;
	}
	if (in_left_files.size() == 0) {
		INFO("No files to read from.");
		return 2;
	}
	//end parsing
	seqan::SeqFileOut out_good(out_good_file.c_str());
	seqan::SeqFileOut out_bad;
	if (vm.count("bad")) {
		open(out_bad, out_bad_file.c_str());
	}
	rm::Read left, right, result;
	rm::SequenceMerger merger;
	int cnt_good = 0, cnt_all = 0;
	for (size_t i = 0; i < in_left_files.size(); ++i) {
		seqan::SeqFileIn in_left(in_left_files[i].c_str());
		seqan::SeqFileIn in_right(in_right_files[i].c_str());
		while (!atEnd(in_left) && !atEnd(in_right)) {
			left.read(in_left);
			right.read(in_right);
			right.reverseRead();
		//	pair <rm::Read, bool> result = merger.Merge<rm::SuffArrayPairer>(left,right, settings, cnt_all);	
		//	pair <rm::Read, bool> result = merger.Merge<rm::SimplePairer>(left,right, settings, cnt_all);	
		//	pair <rm::Read, bool> result = merger.Merge<rm::KHashPairer>(left,right, settings, cnt_all);	
			pair <rm::Read, bool> result = merger.Merge<rm::SmartKHashPairer>(left,right, settings, cnt_all);	
			if (result.second) {
				result.first.write(out_good);
				++cnt_good;
			} else if (vm.count("bad")){
				left.write(out_bad);
				right.write(out_bad);
			}
			++cnt_all;
			if (cnt_all % 1000 == 0) {
				INFO(cnt_all << " reads were processed");
			}
		}
		INFO("Finished reading from " << in_left_files[i] << " and " << in_right_files[i]);
	}
	INFO(cnt_good << " reads from " << cnt_all << " were successfully merged");
	INFO("Merged reads were written to " << out_good_file);
	if (vm.count("bad")) {
		INFO("Non-merged reads were written to " << out_bad_file);
	}
	return 0;
}
