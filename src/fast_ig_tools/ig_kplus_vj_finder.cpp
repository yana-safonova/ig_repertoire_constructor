#include <memory>
#include <cassert>
#include <iostream>
#include <vector>
#include <fstream>
#include <boost/algorithm/string/replace.hpp>
#include <mutex>
#include <stdexcept>

using std::vector;
using std::cout;
using std::cin;
using std::cerr;
using std::endl;

#include <boost/program_options.hpp>
#include "fast_ig_tools.hpp"
using path::make_dirs;

#include <seqan/seq_io.h>
using seqan::Dna5String;
using seqan::CharString;
using seqan::length;
using seqan::SeqFileIn;
using seqan::SeqFileOut;

#include "ig_block_alignment.hpp"
#include "ig_kplus_vj_finder.hpp"

using namespace fast_ig_tools;

using std::string;


struct Ig_KPlus_Finder_Parameters {
    int K = 7; // anchor length
    int word_size_j = 5;
    size_t left_uncovered_limit = 16;
    size_t right_uncovered_limit = 5; // It should be at least 4 (=1 + 3cropped) 1bp trimming is common
    size_t min_v_segment_length = 250;
    size_t min_j_segment_length = 30;
    std::string input_file = "", organism = "human";
    int max_local_deletions = 12;
    int max_local_insertions = 12;
    int max_global_gap = 24;
    int min_k_coverage = 50;
    int min_k_coverage_j = 13;
    int max_candidates = 10;
    int max_candidates_j = 10;
    std::string loci = "all";
    std::string db_directory = "./germline";
    std::string output_dir;
    size_t threads = 4;
    bool verbose = false;
    bool compress = false;
    bool pseudogenes = true;
    bool fix_spaces = true;
    std::string output_filename;
    std::string bad_output_filename;
    std::string add_info_filename;
    std::string discard_info_filename;
    std::string vgenes_filename;
    std::string jgenes_filename;
    const std::string separator_default = "comma";
    std::string separator = "";
    size_t min_len = 300;

    bool fill_left = true;
    bool fill_right = false;
    size_t fix_left = 3;
    size_t fix_right = 0;
    bool crop_left = true;
    bool crop_right = true;

    void parse_cmd_line(int argc, char **argv) {
        namespace po = boost::program_options;
        using std::exit;

        // Declare a group of options that will be
        // allowed only on command line
        po::options_description generic("Generic options");
        generic.add_options()
            ("version,v", "print version string")
            ("help,h", "produce help message")
            ("config-file,c", "name of a file of a configuration")
            ("input-file,i", po::value<std::string>(&input_file)->required(),
             "name of an input file (FASTA|FASTQ)")
            ("output-dir,o", po::value<std::string>(&output_dir)->required(),
             "output directory")
            ;

        auto separator_notifier = [this](const std::string &arg) {
            if (arg == "comma") {
                separator = ",";
            } else if (arg == "semicolon") {
                separator = ";";
            } else if (arg == "tab" || arg == "tabular") {
                separator = "\t";
            } else {
                separator = arg;
            }
        };

        // Declare a group of options that will be
        // allowed both on command line and in
        // config file
        po::options_description config("Configuration");
        config.add_options()
            ("compress,Z", po::value<bool>(&compress)->default_value(compress),
             "compress output FASTA files using zlib")
            ("pseudogenes,P", po::value<bool>(&pseudogenes)->default_value(pseudogenes),
             "use pseudogenes along with normal germline genes")
            ("verbose,V", po::value<bool>(&verbose)->default_value(verbose),
             "produce alignment output for each query")
            ("loci,l", po::value<std::string>(&loci)->default_value(loci),
             "loci: IGH, IGL, IGK, IG (all BCRs), TRA, TRB, TRG, TRD, TR (all TCRs) or all")
            ("db-directory", po::value<std::string>(&db_directory)->default_value(db_directory),
             "directory with germline database")
            ("threads,t", po::value<size_t>(&threads)->default_value(threads),
             "the number of threads")
            ("word-size,k", po::value<int>(&K)->default_value(K),
             "word size for V-genes")
            ("word-size-j", po::value<int>(&word_size_j)->default_value(word_size_j),
             "word size for J-genes")
            ("min-k-coverage,n", po::value<int>(&min_k_coverage)->default_value(min_k_coverage),
             "minimal k+-coverage")
            ("min-k-coverage-j", po::value<int>(&min_k_coverage_j)->default_value(min_k_coverage_j),
             "minimal k+-coverage for J-gene")
            ("max-candidates-v,N", po::value<int>(&max_candidates)->default_value(max_candidates),
             "maximal number of candidates for each query")
            ("max-candidates-j", po::value<int>(&max_candidates_j)->default_value(max_candidates_j),
             "maximal number of J-gene candidates for each query")
            ("organism", po::value<std::string>(&organism)->default_value(organism),
             "organism ('human', 'mouse', 'pig', 'rabbit', 'rat' and 'rhesus_monkey' are supported)")
            ("fix-spaces", po::value<bool>(&fix_spaces)->default_value(fix_spaces),
             "replace spaces in read ids by underline symbol '_'")
            ("separator", po::value<std::string>()->default_value(separator_default)->notifier(separator_notifier),
             "separator for alignment info file: 'comma', 'semicolon', 'tab' (or 'tabular') or custom string")
            ("min-len", po::value<size_t>(&min_len)->default_value(min_len),
             "minimal length of reported sequence")
            ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        po::options_description hidden("Hidden options");
        hidden.add_options()
            ("help-hidden", "show all options, including developers' ones")
            ("left-uncoverage-limit", po::value<size_t>(&left_uncovered_limit)->default_value(left_uncovered_limit),
             "uncoverage limit of left end")
            ("right-uncoverage-limit", po::value<size_t>(&right_uncovered_limit)->default_value(right_uncovered_limit),
             "uncoverage limit of right end")
            ("max-global-gap", po::value<int>(&max_global_gap)->default_value(max_global_gap),
             "maximal allowed size of global gap")
            ("max-local-insertions", po::value<int>(&max_local_insertions)->default_value(max_local_insertions),
             "maximal allowed size of local insertion")
            ("max-local-deletions", po::value<int>(&max_local_deletions)->default_value(max_local_deletions),
             "maximal allowed size of local deletion")
            ("min-vsegment-length", po::value<size_t>(&min_v_segment_length)->default_value(min_v_segment_length),
             "minimal allowed length of V gene segment")
            ("min-jsegment-length", po::value<size_t>(&min_j_segment_length)->default_value(min_j_segment_length),
             "minimal allowed length of J gene segment")
            ("fill-left", po::value<bool>(&fill_left)->default_value(fill_left),
             "fill left cropped positions by germline")
            ("fill-right", po::value<bool>(&fill_right)->default_value(fill_right),
             "fill right cropped positions by germline")
            ("crop-left", po::value<bool>(&crop_left)->default_value(crop_left),
             "crop extra left positions")
            ("crop-right", po::value<bool>(&crop_right)->default_value(crop_right),
             "crop extra right positions")
            ("fix-left", po::value<size_t>(&fix_left)->default_value(fix_left),
             "the number left read positions which will be fixed by germline")
            ("fix-right", po::value<size_t>(&fix_right)->default_value(fix_right),
             "the number right read positions which will be fixed by germline")
            ;

        po::options_description cmdline_options("All command line options");
        cmdline_options.add(generic).add(config).add(hidden);

        po::options_description config_file_options;
        config_file_options.add(config).add(hidden);

        po::options_description visible("Allowed options");
        visible.add(generic).add(config);

        po::positional_options_description p;
        p.add("input-file", 1).add("output-dir", 1);

        po::variables_map vm;
        store(po::command_line_parser(argc, argv).
              options(cmdline_options).positional(p).run(), vm);

        if (vm.count("help-hidden")) {
            cout << cmdline_options << std::endl;
            exit(0);
        }

        if (vm.count("help")) {
            cout << visible << "\n";
            exit(0);
        }

        if (vm.count("version")) {
            cout << "ig_kplus_vjfinder version 0.2" << std::endl;
            exit(0);
        }

        if (vm.count("config-file")) {
            std::string config_file = vm["config-file"].as<std::string>();
            std::ifstream ifs(config_file.c_str());
            if (!ifs) {
                ERROR("Config file " << config_file << " was not found");
                exit(1);
            } else {
                store(parse_config_file(ifs, config_file_options), vm);
                // reparse cmd line again for update config defaults
                store(po::command_line_parser(argc, argv).
                      options(cmdline_options).positional(p).run(), vm);
            }
        }

        try {
            notify(vm);
        } catch (const po::error &e) {
            ERROR("Command-line parser error: " << e.what());
            exit(1);
        } catch (const std::exception &e) {
            ERROR("Unknown exception: " << e.what());
            exit(1);
        }

        // if (separator == "comma") {
        //     separator = ",";
        // } else if (separator == "semicolon") {
        //     separator = ";";
        // } else if (separator == "tab" || separator == "tabular") {
        //     separator = "\t";
        // } else {
        //     // Let it be. Do nothing
        // }
    }

    void print_info() const {
        INFO(bformat("Input FASTQ reads: %s") % input_file);
        INFO(bformat("Output directory: %s") % output_dir);
        INFO(bformat("Organism: %s") % organism);
        INFO(bformat("Locus: %s") % loci);
        INFO("Word size for V-gene: " << K);
        INFO("Word size for J-gene: " << word_size_j);
    }

    void prepare_output() {
        make_dirs(output_dir);
        output_filename = output_dir + "/cleaned_reads.fa";
        bad_output_filename = output_dir + "/filtered_reads.fa";
        add_info_filename = output_dir + "/alignment_info.csv";
        discard_info_filename = output_dir + "/discard_info.txt";

        if (compress) {
            output_filename += ".gz";
            bad_output_filename += ".gz";
        }
    }
};

template<typename T>
size_t replace_spaces(T &s, char target = '_') {
    size_t count = 0;
    for (size_t i = 0, l = length(s); i < l; ++i) {
        if (s[i] == ' ') {
            s[i] = target;
            ++count;
        }
    }

    return count;
}

int main(int argc, char **argv) {
    segfault_handler sh;
    perf_counter pc;
    // create_console_logger("./fast_ig_tools.log");
    create_console_logger("");

    INFO("Command line: " << join_cmd_line(argc, argv));

    Ig_KPlus_Finder_Parameters param;
    param.parse_cmd_line(argc, argv);
    param.prepare_output();

    const VJAligner db(param);

    param.print_info();
    INFO("V gene germline database size: " << db.all_loci_database.v_reads.size());
    INFO("J gene germline database size: " << db.all_loci_database.j_reads.size());

    seqan::SeqFileIn seqFileIn_reads(param.input_file.c_str());

    std::mutex stdout_mtx;

    vector<CharString> read_ids;
    vector<Dna5String> reads;
    readRecords(read_ids, reads, seqFileIn_reads);
    INFO(reads.size() << " reads were extracted from " << param.input_file);

    // Fix spaces if asked
    if (param.fix_spaces) {
        for (auto &id : read_ids) {
            replace_spaces(id);
        }
    }

    vector<int> output_isok(reads.size());  // Do not use vector<bool> here due to it is not thread-safe
    std::vector<std::string> add_info_strings(reads.size());
    std::string output_pat = "%s, %d, %d, %1.2f, %s, %d, %d, %1.2f, %s";
    boost::replace_all(output_pat, ", ", param.separator);

    omp_set_num_threads(param.threads);
    INFO(bformat("Alignment reads using %d threads starts") % param.threads);

    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 8))
    for (size_t j = 0; j < reads.size(); ++j) {
        const CharString &read_id = read_ids[j];
        output_isok[j] = false;

        if (length(reads[j]) < param.min_len) {
            // Discard read
            bformat bf("Read too short: %d");
            bf % length(reads[j]);
            add_info_strings[j] = bf.str();
            continue;
        }

        const auto RESULT = db.Query(reads[j], true, true, param);
        if (param.verbose) {
            std::lock_guard<std::mutex> lck(stdout_mtx); //TODO Use OMP critical section

            cout << "Query: " << read_id << endl;
            cout << "Strand: " << RESULT.Strand() << endl;
            cout << "V genes:" << endl;

            if (RESULT.VHitsSize()) {
                for (size_t i = 0; i < RESULT.VHitsSize(); ++i) {
                    cout << bformat("k+-score %3%; %1%:%2%; V-gene: %4%\n")
                        % (RESULT.VStart(i) + 1 ) % RESULT.VEnd(i)
                        % RESULT.VScore(i)        % RESULT.VId(i);

                    cout << RESULT.VAlignmentSeqAn(i);
                }
            } else {
                cout << "V genes not found" << endl;
            }

            if (RESULT.JHitsSize()) {
                for (size_t i = 0; i < RESULT.JHitsSize(); ++i) {
                    cout << bformat("k+-score %3%; %1%:%2%; J-gene: %4%\n")
                        % (RESULT.JStart(i) + 1 )  % RESULT.JEnd(i)
                        % RESULT.JScore(i)         % RESULT.JId(i);

                    cout << RESULT.JAlignmentSeqAn(i);
                }
            } else {
                cout << "J genes not found" << endl;
            }
        }

        if (!RESULT.VHitsSize()) {
            // Discard
            add_info_strings[j] = "No V genes found";
            continue;
        }

        if (!RESULT.JHitsSize()) {
            // Discard
            add_info_strings[j] = "No J genes found";
            continue;
        }

        if (RESULT.RightUncovered() > param.right_uncovered_limit) {
            // Discard read
            bformat bf("Right cropped: %d");
            bf % RESULT.RightUncovered();
            add_info_strings[j] = bf.str();
            continue;
        }

        if (RESULT.LeftUncovered() > param.left_uncovered_limit) {
            // Discard read
            bformat bf("Left cropped: %d");
            bf % RESULT.LeftUncovered();
            add_info_strings[j] = bf.str();
            continue;
        }

        if (RESULT.VSegmentLength() < param.min_v_segment_length) {
            // Discard read
            bformat bf("V segment is too short: %d");
            bf % RESULT.VSegmentLength();
            add_info_strings[j] = bf.str();
            continue;
        }

        if (RESULT.JSegmentLength() < param.min_j_segment_length) {
            // Discard read
            bformat bf("J segment is too short: %d");
            bf % RESULT.JSegmentLength();
            add_info_strings[j] = bf.str();
            continue;
        }

        if (RESULT.FinalLength() < param.min_len) {
            // Discard read
            bformat bf("Final sequence is too short: %d");
            bf % RESULT.FinalLength();
            add_info_strings[j] = bf.str();
            continue;
        }

        bformat bf(output_pat);
        bf  % read_id
            % (RESULT.VStart() + 1) % RESULT.VEnd()
            % RESULT.VScore()       % RESULT.VId()
            % (RESULT.JStart() + 1) % RESULT.JEnd()
            % RESULT.JScore()       % RESULT.JId();

        add_info_strings[j] = bf.str();

        reads[j] = RESULT.FixCropFill(param);
        output_isok[j] = true;
    }

    INFO("Saving results");
    seqan::SeqFileOut cropped_reads_seqFile(param.output_filename.c_str());
    seqan::SeqFileOut bad_reads_seqFile(param.bad_output_filename.c_str());
    std::ofstream add_info(param.add_info_filename.c_str());
    std::ofstream discard_info(param.discard_info_filename.c_str());
    std::string pat = "%s, %s, %s, %s, %s, %s, %s, %s, %s\n";
    boost::replace_all(pat, ", ", param.separator);
    add_info << bformat(pat)
        % "id"
        % "Vstart" % "Vend"
        % "Vscore" % "Vid"
        % "Jstart" % "Jend"
        % "Jscore" % "Jid";

    size_t good_reads = 0;
    for (size_t j = 0; j < reads.size(); ++j) {
        if (output_isok[j]) {
            seqan::writeRecord(cropped_reads_seqFile, read_ids[j], reads[j]);
            add_info << add_info_strings[j] << "\n";
            ++good_reads;
        } else {
            seqan::writeRecord(bad_reads_seqFile, read_ids[j], reads[j]);
            discard_info << read_ids[j] << "\n";
            discard_info << add_info_strings[j] << "\n";
        }
    }

    INFO("Alignment reads finished");
    size_t num_good_reads = good_reads;
    size_t num_bad_reads = reads.size() - good_reads;
    INFO(num_good_reads << " Ig-Seq reads were written to " << param.output_filename);
    INFO(num_bad_reads << " junk (not Ig-Seq) reads were written to " << param.bad_output_filename);
    INFO("Running time: " << running_time_format(pc));
    return 0;
}

// vim: ts=4:sw=4
