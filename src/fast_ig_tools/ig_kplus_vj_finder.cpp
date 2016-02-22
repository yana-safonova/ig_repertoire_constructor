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
    int left_uncoverage_limit = 16;
    int right_uncoverage_limit = 5; // It should be at least 4 (=1 + 3cropped) 1bp trimming is common
    int min_vsegment_length = 250;
    int min_jsegment_length = 30;
    std::string input_file = "", organism = "human";
    int max_local_deletions = 12;
    int max_local_insertions = 12;
    int max_global_gap = 24;
    int min_k_coverage = 50;
    int min_k_coverage_j = 13;
    int max_candidates = 3;
    int max_candidates_j = 3;
    std::string loci = "all";
    std::string db_directory = "./germline";
    std::string output_dir;
    int threads = 4;
    bool silent = true;
    bool fill_prefix_by_germline = true;
    bool compress = false;
    bool pseudogenes = true;
    bool fix_spaces = true;
    std::string config_file = "";
    std::string output_filename;
    std::string bad_output_filename;
    std::string add_info_filename;
    std::string vgenes_filename;
    std::string jgenes_filename;
    int left_fill_germline = 3;
    std::string separator = "comma";
    size_t min_len = 300;

    std::unique_ptr<GermlineDB> db;

    Ig_KPlus_Finder_Parameters(const Ig_KPlus_Finder_Parameters &) = delete;
    Ig_KPlus_Finder_Parameters(Ig_KPlus_Finder_Parameters &&) = delete;
    Ig_KPlus_Finder_Parameters& operator=(const Ig_KPlus_Finder_Parameters &) = delete;
    Ig_KPlus_Finder_Parameters& operator=(Ig_KPlus_Finder_Parameters &&) = delete;

    explicit Ig_KPlus_Finder_Parameters(int argc, char **argv) {
        namespace po = boost::program_options;
        using std::exit;

        // Declare a group of options that will be
        // allowed only on command line
        po::options_description generic("Generic options");
        generic.add_options()
            ("version,v", "print version string")
            ("help,h", "produce help message")
            ("config,c", po::value<std::string>(&config_file)->default_value(config_file),
             "name of a file of a configuration")
            ("input-file,i", po::value<std::string>(&input_file),
             "name of an input file (FASTA|FASTQ)")
            ;

        // Declare a group of options that will be
        // allowed both on command line and in
        // config file
        po::options_description config("Configuration");
        config.add_options()
            ("output,o", po::value<std::string>(&output_dir),
             "output directory")
            ("compress,Z", "compress FASTA files using Zlib")
            ("no-compress", "don't compress FASTA files (default)")
            ("pseudogenes,P", "use pseudogenes along with normal genes (default)")
            ("no-pseudogenes", "don't use pseudogenes along with normal genes")
            ("silent,S", "suppress output for each query (default)")
            ("no-silent,V", "produce info output for each query")
            ("loci,l", po::value<std::string>(&loci)->default_value(loci),
             "loci: IGH, IGL, IGK, IG (all BCRs), TRA, TRB, TRG, TRD, TR (all TCRs) or all")
            ("db-directory", po::value<std::string>(&db_directory)->default_value(db_directory),
             "directory with germline database")
            ("threads,t", po::value<int>(&threads)->default_value(threads),
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
             "organism ('human', 'mouse', 'pig', 'rabbit', 'rat' and 'rhesus_monkey' are supported for this moment)")
            ("fill-prefix-by-germline",
             "fill truncated V-gene prefix by germline content")
            ("no-fill-prefix-by-germline (default)",
             "fill truncated V-gene prefix by 'N'-s")
            ("fix-spaces",
             "replace spaces in read ids by underline symbol '_' (default)")
            ("no-fix-spaces",
             "save spaces in read ids, do not replace them")
            ("separator", po::value<std::string>(&separator)->default_value(separator),
             "separator for alignment info file: 'comma' (default), 'semicolon', 'tab' (or 'tabular') or custom string)")
            ("min-len", po::value<size_t>(&min_len)->default_value(min_len),
             "minimal length of reported sequence")
            ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        po::options_description hidden("Hidden options");
        hidden.add_options()
            ("help-hidden", "show all options, including developers' ones")
            ("left-uncoverage-limit", po::value<int>(&left_uncoverage_limit)->default_value(left_uncoverage_limit),
             "uncoverage limit of left end")
            ("right-uncoverage-limit", po::value<int>(&right_uncoverage_limit)->default_value(right_uncoverage_limit),
             "uncoverage limit of right end")
            ("max-global-gap", po::value<int>(&max_global_gap)->default_value(max_global_gap),
             "maximal allowed size of global gap")
            ("max-local-insertions", po::value<int>(&max_local_insertions)->default_value(max_local_insertions),
             "maximal allowed size of local insertion")
            ("max-local-deletions", po::value<int>(&max_local_deletions)->default_value(max_local_deletions),
             "maximal allowed size of local deletion")
            ("left-fill-germline", po::value<int>(&left_fill_germline)->default_value(left_fill_germline),
             "the number left positions which will be filled by germline")
            ("min-vsegment-length", po::value<int>(&min_vsegment_length)->default_value(min_vsegment_length),
             "minimal allowed length of V gene segment")
            ("min-jsegment-length", po::value<int>(&min_jsegment_length)->default_value(min_jsegment_length),
             "minimal allowed length of J gene segment")
            ;

        po::options_description cmdline_options("All command line options");
        cmdline_options.add(generic).add(config).add(hidden);

        po::options_description config_file_options;
        config_file_options.add(config).add(hidden);

        po::options_description visible("Allowed options");
        visible.add(generic).add(config);

        po::positional_options_description p;
        p.add("input-file", -1);

        po::variables_map vm;
        store(po::command_line_parser(argc, argv).
              options(cmdline_options).positional(p).run(), vm);
        notify(vm);

        if (config_file != "") {
            std::ifstream ifs(config_file.c_str());
            if (!ifs) {
                ERROR("Config file " << config_file << " was not found");
                exit(1);
            } else {
                store(parse_config_file(ifs, config_file_options), vm);
                // reparse cmd line again for update config defaults
                store(po::command_line_parser(argc, argv).
                      options(cmdline_options).positional(p).run(), vm);
                notify(vm);
            }
        }

        if (vm.count("help-hidden")) {
            cout << cmdline_options << std::endl;
            exit(0);
        }

        if (vm.count("help") || !vm.count("input-file")) { // TODO Process required arguments by the proper way
            cout << visible << "\n";
            exit(0);
        }

        if (vm.count("version")) {
            cout << "<Some cool name> version 0.1" << vm.count("version") << std::endl;
            exit(0);
        }

        if (vm.count("silent")) {
            silent = true;
        } else if (vm.count("no-silent")) {
            silent = false;
        }

        if (vm.count("compress")) {
            compress = true;
        } else if (vm.count("no-compress")) {
            compress = false;
        }

        if (vm.count("pseudogenes")) {
            pseudogenes = true;
        } else if (vm.count("no-pseudogenes")) {
            pseudogenes = false;
        }

        if (vm.count("no-fill-prefix-by-germline")) {
            fill_prefix_by_germline = false;
        }

        if (vm.count("fill-prefix-by-germline")) {
            fill_prefix_by_germline = true;
        }

        if (vm.count("no-fix-spaces")) {
            fix_spaces = false;
        }

        if (vm.count("fix-spaces")) {
            fix_spaces = true;
        }

        if (separator == "comma") {
            separator = ",";
        } else if (separator == "semicolon") {
            separator = ";";
        } else if (separator == "tab" || separator == "tabular") {
            separator = "\t";
        } else {
            // Let it be. Do nothing
        }

        INFO(bformat("Input FASTQ reads: %s") % input_file);
        INFO(bformat("Output directory: %s") % output_dir);
        INFO(bformat("Organism: %s") % organism);
        INFO(bformat("Locus: %s") % loci);
        INFO("Word size for V-gene: " << K);
        INFO("Word size for J-gene: " << word_size_j);

        prepare_output();

        db.reset(new GermlineDB(*this));

        INFO("V-gene germline database size: " << db->all_v_reads.size());
        INFO("J-gene germline database size: " << db->all_j_reads.size());
    }

private:
    void prepare_output() {
        make_dirs(output_dir);
        output_filename = output_dir + "/cleaned_reads.fa";
        bad_output_filename = output_dir + "/filtered_reads.fa";
        add_info_filename = output_dir + "/alignment_info.csv";

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

    Ig_KPlus_Finder_Parameters param(argc, argv);

    const auto &v_reads = param.db->all_v_reads;
    const auto &j_reads = param.db->all_j_reads;
    const auto &v_ids = param.db->all_v_ids;
    const auto &j_ids = param.db->all_j_ids;

    seqan::SeqFileIn seqFileIn_reads(param.input_file.c_str());

    std::mutex stdout_mtx;
    const BlockAligner index(v_reads, param.K, param.max_global_gap, param.left_uncoverage_limit,
                             1005000,
                             param.max_local_insertions, param.max_local_deletions, param.min_k_coverage);
    const BlockAligner j_index(j_reads, param.word_size_j,
                               param.max_global_gap, 100000, 10000,
                               param.max_local_insertions, param.max_local_deletions, param.min_k_coverage_j);

    vector<CharString> read_ids;
    vector<Dna5String> reads;
    readRecords(read_ids, reads, seqFileIn_reads);
    INFO(reads.size() << " reads were extracted from " << param.input_file);

    // Fix spaces if asked
    for (auto &id : read_ids) {
        replace_spaces(id);
    }

    vector<int> output_isok(reads.size());  // Do not use vector<bool> here due to it is not thread-safe
    std::vector<std::string> add_info_strings(reads.size());
    std::string output_pat = "%s, %d, %d, %1.2f, %s, %d, %d, %1.2f, %s";
    boost::replace_all(output_pat, ", ", param.separator);

    omp_set_num_threads(param.threads);
    INFO(bformat("Alignment reads using %d threads starts") % param.threads);

    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 8))
    for (size_t j = 0; j < reads.size(); ++j) {
        CharString read_id = read_ids[j];

        // temporary solution - we crop three nucleotides at the end of string
        if (length(reads[j]) < param.min_len) {
            // Discard so short read
            output_isok[j] = false;
            continue;
        }

        auto RESULT = param.db->query(reads[j], true, true, param);

        Dna5String read = reads[j];
        Dna5String read_rc = read;
        reverseComplement(read_rc);

        auto result_pstrand = index.query(read, param.max_candidates);
        auto result_nstrand = index.query(read_rc, param.max_candidates);

        int pscore = (result_pstrand.size() > 0) ? result_pstrand[0].kp_coverage : 0;
        int nscore = (result_nstrand.size() > 0) ? result_nstrand[0].kp_coverage : 0;

        int strand = (pscore >= nscore) ? 1 : -1;
        const auto &stranded_read = (strand == 1) ? read : read_rc;
        const auto &result = (strand == 1) ? result_pstrand : result_nstrand;

        bool aligned = false;
        if (!result.empty()) { // If we found at least one alignment
            for (const auto &align : result) {
                if (-align.start > param.left_uncoverage_limit) {
                    // Discard read
                    break;
                }

                if (static_cast<int>(align.last_match_read_pos()) - align.start < param.min_vsegment_length) { // TODO check +-1
                    // Discard read
                    break;
                }

                auto last_match = align.path[align.path.size() - 1];
                int end_of_v = last_match.read_pos + last_match.length;
                auto suff = suffix(stranded_read, end_of_v);
                auto jresult = j_index.query(suff, param.max_candidates_j);

                if (!param.silent) {
                    std::lock_guard<std::mutex> lck(stdout_mtx);

                    cout << "Query: " << read_id << endl;
                    cout << bformat("Strand: %s (best V-gene k+-score: +%d/-%d)")
                        % ((strand == 1) ? "+" : "-") % pscore % nscore <<  endl; // '?' has less priority

                    cout << "V-genes:" << endl;

                    cout << bformat("k+-score %3%; %1%:%2%; V-gene: %4%")
                        % (align.start+1)   % end_of_v
                        % align.kp_coverage % v_ids[align.needle_index];

                    cout << " " << align.path.visualize_matches(length(v_reads[align.needle_index]), length(stranded_read)) << endl;
                    if (!jresult.empty()) {
                        cout << "\tJ-genes:" << endl;
                        for (const auto &jalign : jresult) {
                            cout << bformat("\tk+-coverage %3%; %1%:%2%; J-gene: %4%")
                                % (jalign.start+1 + end_of_v) % (jalign.finish + end_of_v)
                                % jalign.kp_coverage          % j_ids[jalign.needle_index];

                            cout << " " << jalign.path.visualize_matches(length(j_reads[jalign.needle_index]), length(suff)) << endl;
                        }
                    } else {
                        cout << "J-genes not found" << endl;
                    }
                    cout << endl;
                }

                if (!jresult.empty() && !aligned) { // Save as <<good>> read
                    Dna5String cropped_read; // Crop to head of V-gene
                    if (align.start >= 0) {
                        cropped_read = suffix(stranded_read, align.start);
                    } else {
                        if (param.fill_prefix_by_germline) {
                            cropped_read = prefix(v_reads[align.needle_index], std::abs(align.start)); // Fill by germline
                        } else {
                            cropped_read = std::string(std::abs(align.start), 'N');
                        }
                        cropped_read += stranded_read;
                    }

                    if (param.left_fill_germline > 0) {
                        for (int i = 0; i < param.left_fill_germline; ++i) {
                            cropped_read[i] = v_reads[align.needle_index][i]; // Fill by germline
                        }
                    }

                    const auto &jalign = *jresult.cbegin();

                    if (static_cast<int>(jalign.finish - jalign.first_match_read_pos()) < param.min_jsegment_length) { // TODO check +-1
                        // Discard read
                        break;
                    }

                    if (jalign.finish - int(length(suff)) > param.right_uncoverage_limit) {
                        // Discard read
                        // INFO(suff);
                        // INFO("Read discarded by J gene threshold " << jalign.finish << " " << length(suff));
                        break;
                    }

                    {
                        const auto &first_jalign = jalign.path[0];

                        bformat bf(output_pat);
                        bf % read_id
                           % (align.start+1)             % end_of_v
                           % align.score                 % toCString(v_ids[align.needle_index])
                           % (first_jalign.read_pos + 1 + end_of_v) % (jalign.finish + end_of_v)
                           % jalign.score                % toCString(j_ids[jalign.needle_index]);

                        add_info_strings[j] = bf.str();

                        if (length(cropped_read) >= param.min_len) {
                            reads[j] = cropped_read;
                            output_isok[j] = true;
                        } else {
                            // Read is too short
                            output_isok[j] = false;
                        }
                    }

                    aligned = true;
                }
            }
        } else if (!param.silent) {
            std::lock_guard<std::mutex> lck(stdout_mtx);

            cout << "Query: " << read_id << endl;
            cout << bformat("Strand: %s (best V-gene k+-score: +%d/-%d)")
                % ((strand == 1) ? "+" : "-") % pscore % nscore <<  endl; // '?' has less priority
            cout << "V-genes not found" << endl;
            cout << endl;
        }

        if (!aligned) { // Save as bad read
            // output_reads[j] = read;
            // Do nothing, save origonal read
            output_isok[j] = false;
        }
    }

    INFO("Saving results");
    seqan::SeqFileOut cropped_reads_seqFile(param.output_filename.c_str());
    seqan::SeqFileOut bad_reads_seqFile(param.bad_output_filename.c_str());
    std::ofstream add_info(param.add_info_filename.c_str());
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
