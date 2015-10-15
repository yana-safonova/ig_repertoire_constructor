#include <seqan/seq_io.h>
#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <memory>
#include <unordered_map>
#include <fstream>
#include <boost/format.hpp>
#include <mutex>
#include <chrono>

using seqan::Dna5String;
using seqan::CharString;
using seqan::length;
using seqan::SeqFileIn;
using seqan::SeqFileOut;

using std::vector;
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::map;
using std::make_pair;
using bformat = boost::format;

#include <boost/program_options.hpp>
#include "fast_ig_tools.hpp"
using path::make_dirs;


namespace fast_ig_tools {

// input vector(pair(needle_pos, read_pos))
// output vector(tuple(needle_pos, read_pos, length))
struct Match {
    int needle_pos;
    int read_pos;
    size_t length;
};


struct KmerMatch {
    int needle_pos;
    int read_pos;
};


int path_coverage_length(const std::vector<Match> &path) {
    int result = 0;

    for (const auto &match : path) {
        result += match.length;
    }

    return result;
}


std::vector<Match> combine_sequential_kmer_matches(const std::vector<KmerMatch> &kmatches,
                                                   size_t K) {
    // Convert to (shift, read_pos)
    std::vector<std::pair<int, int>> matches;
    matches.reserve(kmatches.size());
    for (auto &kmatch : kmatches) {
        matches.push_back( { kmatch.read_pos - kmatch.needle_pos, kmatch.read_pos } );
    }

    std::sort(matches.begin(), matches.end());

    std::vector<Match> res;
    res.reserve(matches.size());

    if (matches.size() == 0) {
        return res;
    }

    Match cur = { matches[0].second - matches[0].first, matches[0].second, K }; // start first match
    for (size_t i = 1; i < matches.size(); ++i) {
        if (matches[i].first == matches[i-1].first && matches[i].second == matches[i-1].second + 1) { // extend current match
            cur.length += 1;
        } else { // save match and start new one
            res.push_back(cur);
            cur = { matches[i].second - matches[i].first, matches[i].second, K };
        }
    }
    res.push_back(cur); // save last match

    return res;
}


std::string visualize_matches(const std::vector<Match> &matches,
                              int needle_length, int read_length) {
    // Draw fancy alignment
    // (read/needle)
    std::stringstream ss;

    assert(std::is_sorted(matches.cbegin(), matches.cend(),
                          [](const Match &a, const Match &b) -> bool { return a.needle_pos < b.needle_pos; }));
    assert(std::is_sorted(matches.cbegin(), matches.cend(),
                          [](const Match &a, const Match &b) -> bool { return a.read_pos < b.read_pos; }));

    ss << bformat("{%d}") % std::max(matches[0].needle_pos - matches[0].read_pos, 0);
    ss << bformat("(%d)") % std::min(matches[0].needle_pos, matches[0].read_pos);
    for (size_t i = 0; i < matches.size() - 1; ++i) {
        int read_gap = matches[i+1].read_pos - matches[i].read_pos - matches[i].length;
        int needle_gap = matches[i+1].needle_pos - matches[i].needle_pos - matches[i].length;
        unsigned current_match_len = matches[i].length;

        if (needle_gap >=0 || read_gap >= 0) {
            if (needle_gap != read_gap) {
                ss << bformat("%d(%d%+d)") % current_match_len % needle_gap % (read_gap - needle_gap);
            } else {
                ss << bformat("%2%(%1%)") % read_gap % current_match_len;
            }
        }
    }

    const auto &last_match = matches[matches.size() - 1];
    ss << bformat("%1%") % last_match.length;
    ss << bformat("(%d)") % std::min(needle_length - last_match.needle_pos, read_length - last_match.read_pos);
    ss << bformat("{%d}") % std::max((needle_length - last_match.needle_pos) - (read_length - last_match.read_pos), 0);

    return ss.str();
}


class KmerIndex {
    struct PositionInDB {
        size_t needle_index;
        size_t position;
    };

    struct Alignment {
        int kp_coverage;
        std::vector<Match> path;
        int start, finish;
        size_t needle_index;
        int overlap_length;
        double score;

        bool operator< (const Alignment& b) const {
            return this->kp_coverage < b.kp_coverage;
        }
    };

public:
    KmerIndex(const vector<Dna5String> &queries,
              size_t K,
              int max_global_gap, int left_uncoverage_limit, int right_uncoverage_limit,
              int max_local_insertions,
              int max_local_deletions,
              int min_k_coverage) : max_local_insertions{max_local_insertions},
                                    max_local_deletions{max_local_deletions},
                                    min_k_coverage{min_k_coverage},
                                    K{K},
                                    left_uncoverage_limit{left_uncoverage_limit},
                                    right_uncoverage_limit{right_uncoverage_limit},
                                    max_global_gap{max_global_gap} {
            this->queries.resize(queries.size());
            for (size_t i = 0; i < queries.size(); ++i) {
                this->queries[i] = queries[i];
            }

            for (size_t j = 0; j < length(queries); ++j) {
                for (size_t start = 0; start + K <= length(queries[j]); ++start) {
                    kmer2needle[infix(queries[j], start, start + K)].push_back( { j, start } );
                }
            }
        }

    std::vector<Alignment> query(Dna5String read, size_t limit) const {
        auto result = query_unordered(read);

        using ctuple_type = decltype(*result.cbegin());

        auto score_function = [](ctuple_type a) { return a.kp_coverage; };
        auto comp = [&score_function](ctuple_type a, ctuple_type b) -> bool { return score_function(a) > score_function(b); };

        // Return top <limit>
        std::nth_element(result.begin(), result.begin() + limit, result.end(), comp);
        result.resize(std::min(result.size(), limit));
        std::sort(result.begin(), result.end(), comp);

        return result;
    }

private:
    std::unordered_map<size_t, std::vector<KmerMatch>> needle2matches(Dna5String read) const {
        std::unordered_map<size_t, std::vector<KmerMatch>> result;

        if (length(read) < K) {
            // Return empty result
            return result;
        }

        for (size_t j = 0; j < length(read) - K + 1; ++j) {
            auto kmer = infixWithLength(read, j, K);
            auto it = kmer2needle.find(kmer);

            if (it == kmer2needle.cend()) {
                continue;
            }

            for (const auto &p : it->second) {
                size_t needle_index = p.needle_index;
                size_t kmer_pos_in_read = j;
                size_t kmer_pos_in_needle = p.position;
                int shift = static_cast<int>(kmer_pos_in_read) - static_cast<int>(kmer_pos_in_needle);

                // We make these limits less strict because of possibility of indels
                int shift_min = -left_uncoverage_limit - max_global_gap;
                int shift_max = static_cast<int>(length(read)) - static_cast<int>(length(queries[needle_index])) + right_uncoverage_limit + max_global_gap;

                if (shift >= shift_min && shift <= shift_max) {
                    result[needle_index].push_back( { static_cast<int>(kmer_pos_in_needle), static_cast<int>(kmer_pos_in_read) } );
                }
            }
        }

        return result;
    }

    std::vector<Alignment> query_unordered(Dna5String read) const {
        std::vector<Alignment> result;
        auto n2m = this->needle2matches(read);

        if (!n2m.empty()) { // if we found at least one proper k-mer match
            for (const auto &p : n2m) {
                const auto &needle_index = p.first;
                const auto &matches = p.second;

                std::vector<Match> combined = combine_sequential_kmer_matches(matches, K);

                assert(combined.size() > 0);

                std::sort(combined.begin(), combined.end(),
                          [](const Match &a, const Match &b) -> bool { return a.needle_pos < b.needle_pos; });

                std::vector<double> values(combined.size(), 0.);
                std::vector<size_t> next(combined.size());
                std::iota(next.begin(), next.end(), 0);

                // auto has_edge_old = [](const Match &a, const Match &b) -> bool {
                //   return std::min(b.needle_pos - a.needle_pos, b.read_pos - a.read_pos) >= int(a.length); // signed/unsigned comparisson
                // };
                //
                auto has_edge = [&](const Match &a, const Match &b) -> bool {
                    // int sha = a.read_pos - a.needle_pos;
                    // int shb = b.read_pos - b.needle_pos;

                    int read_gap = b.read_pos - a.read_pos;
                    int needle_gap = b.needle_pos - a.needle_pos;
                    int insert_size = read_gap - needle_gap;

                    if (insert_size > max_local_insertions || -insert_size > max_local_deletions) return false;
                    return std::min(b.needle_pos - a.needle_pos, b.read_pos - a.read_pos) >= int(a.length); // signed/unsigned comparisson
                };

                for (int i = combined.size() - 1; i >= 0; --i) {
                    for (size_t j = i + 1; j < combined.size(); ++j) {
                        if (has_edge(combined[i], combined[j]) && (values[j] > values[i])) {
                            values[i] = values[j];
                            next[i] = j;
                        }
                    }

                    values[i] += combined[i].length;
                }

                std::vector<Match> path;
                path.reserve(combined.size());

                size_t maxi = std::max_element(values.cbegin(), values.cend()) - values.cbegin();

                while (true) {
                    path.push_back(combined[maxi]);
                    if (next[maxi] == maxi) {
                        break;
                    } else {
                        maxi = next[maxi];
                    }
                }

                int coverage_length = path_coverage_length(path);

                assert(std::is_sorted(path.cbegin(), path.cend(), has_edge));

                // Just use the most left and most right matches
                int left_shift = path[0].read_pos - path[0].needle_pos;
                int right_shift = path[path.size() - 1].read_pos - path[path.size() - 1].needle_pos;

                if (std::abs(left_shift - right_shift) > max_global_gap) {
                    // Omit such match
                    continue;
                }

                if (coverage_length < min_k_coverage) {
                    // Omit such match
                    continue;
                }

                int start = left_shift;
                int finish = right_shift + int(length(queries[needle_index]));

                int shift_min = -left_uncoverage_limit;
                int shift_max = int(length(read)) - int(length(queries[needle_index])) + right_uncoverage_limit;

                if (left_shift < shift_min || right_shift > shift_max) {
                    // Omit candidates with unproper final shift
                    // Maybe we should make these limits less strict because of possibility of indels on edges?
                    continue;
                }

                int over_start = std::max(0, start);
                int over_finish = std::min(right_shift + length(queries[needle_index]), length(read));
                int read_overlap_length = over_finish - over_start; // read overlap
                int needle_overlap_length = read_overlap_length + left_shift - right_shift;

                Alignment align;
                align.kp_coverage = coverage_length;
                align.score = static_cast<double>(coverage_length) / static_cast<double>(length(queries[needle_index]));
                align.path = std::move(path);
                align.start = start;
                align.finish = finish;
                align.needle_index = needle_index;
                align.overlap_length = needle_overlap_length;

                result.push_back(std::move(align));
            }
        }

        return result;
    }

    vector<Dna5String> queries;
    int max_local_insertions;
    int max_local_deletions;
    int min_k_coverage;
    size_t K;
    int left_uncoverage_limit, right_uncoverage_limit;
    int max_global_gap;
    map<Dna5String, vector<PositionInDB>> kmer2needle;
};

} // namespace fast_ig_tools

using namespace fast_ig_tools;

using std::string;


struct Ig_KPlus_Finder_Parameters {
    int K = 7; // anchor length
    int word_size_j = 5;
    int left_uncoverage_limit = 16;
    int right_uncoverage_limit = 0; // We have to cover (D)J region too. Maybe it should be even negative
    std::string input_file = "", organism = "human";
    int max_local_deletions = 12;
    int max_local_insertions = 12;
    int max_global_gap = 24;
    int min_k_coverage = 50;
    int min_k_coverage_j = 9;
    int max_candidates = 3;
    int max_candidates_j = 3;
    std::string chain = "all";
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
    size_t num_cropped_nucls = 3;
    std::vector<CharString> v_ids;
    std::vector<Dna5String> v_reads;
    std::vector<CharString> j_ids;
    std::vector<Dna5String> j_reads;

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
            ("chain,C", po::value<std::string>(&chain)->default_value(chain),
             "IG chain: 'all' (default), 'heavy', 'lambda', 'kappa' or 'light' ('lambda' + 'kappa')")
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
             "organism ('human' and 'mouse' are supported for this moment)")
            ("fill-prefix-by-germline",
             "fill truncated V-gene prefix by germline content")
            ("no-fill-prefix-by-germline (default)",
             "fill truncated V-gene prefix by 'N'-s")
            ("fix-spaces",
             "replace spaces in read ids by underline symbol '_' (default)")
            ("no-fix-spaces",
             "save spaces in read ids, do not replace them")
            ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        po::options_description hidden("Hidden options");
        hidden.add_options()
            ("help-hidden", "show all options, including developers' ones")
            ("left-uncoverage-limit", po::value<int>(&left_uncoverage_limit)->default_value(left_uncoverage_limit),
             "uncoverage limit of left end")
            ("max-global-gap", po::value<int>(&max_global_gap)->default_value(max_global_gap),
             "maximal allowed size of global gap")
            ("max-local-insertions", po::value<int>(&max_local_insertions)->default_value(max_local_insertions),
             "maximal allowed size of local insertion")
            ("max-local-deletions", po::value<int>(&max_local_deletions)->default_value(max_local_deletions),
             "maximal allowed size of local deletion")
            ("left-fill-germline", po::value<int>(&left_fill_germline)->default_value(left_fill_germline),
             "the number left positions which will be filled by germline")
            ("right-cropped", po::value<size_t>(&num_cropped_nucls)->default_value(num_cropped_nucls),
             "the number of right positions which will be cropped")
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

        INFO(bformat("Input FASTQ reads: %s") % input_file);
        INFO(bformat("Output directory: %s") % output_dir);
        INFO("Short k-mer size: " << K);

        prepare_output();

        read_genes();
        if (pseudogenes) {
            // Read pseudogenes
            read_genes(true);
        }
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

    std::string gene_file_name(const std::string &chain_letter,
                               const std::string &gene,
                               bool pseudo = false) {
        return db_directory + "/" + organism + "/IG" + chain_letter + gene + (pseudo ? "-pseudo" : "") + ".fa";
    }

    std::vector<std::string> chain_letters() {
        if (chain == "heavy") {
            return { "H" };
        } else if (chain == "lambda") {
            return { "L" };
        } else if (chain == "kappa") {
            return { "K" };
        } else if (chain == "light") {
            return { "K", "L" };
        } else if (chain == "all") {
            return { "K", "L", "H" };
        } else {
            throw "Unknown chain type";
            return {};
        }
    }

    void read_genes(bool pseudo = false) {
        for (const auto &letter : chain_letters()) {
            std::string v_file = gene_file_name(letter, "V", pseudo);
            read_gene(v_file, v_ids, v_reads);

            std::string j_file = gene_file_name(letter, "J", pseudo);
            read_gene(j_file, j_ids, j_reads);
        }
    }

    static void read_gene(const std::string &filename,
                          vector<CharString> &ids, vector<Dna5String> &reads) {
        vector<CharString> _ids;
        vector<Dna5String> _reads;

        SeqFileIn seqFileIn(filename.c_str());
        readRecords(_ids, _reads, seqFileIn);

        ids.reserve(ids.size() + _ids.size());
        for (const auto &id : _ids) {
            ids.push_back(id);
        }

        reads.reserve(reads.size() + _reads.size());
        for (const auto &read : _reads) {
            reads.push_back(read);
        }
    }
};

template<typename T>
size_t replace_spaces(T &s, char target = '_') {
    size_t count = 0;
    for(size_t i = 0, l = length(s); i < l; ++i) {
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

    const auto &v_reads = param.v_reads;
    const auto &j_reads = param.j_reads;
    const auto &v_ids = param.v_ids;
    const auto &j_ids = param.j_ids;

    seqan::SeqFileIn seqFileIn_reads(param.input_file.c_str());

    std::mutex stdout_mtx;
    const KmerIndex index(v_reads, param.K, param.max_global_gap, param.left_uncoverage_limit, param.right_uncoverage_limit,
                          param.max_local_insertions, param.max_local_deletions, param.min_k_coverage);
    const KmerIndex j_index(j_reads, param.word_size_j,
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

    omp_set_num_threads(param.threads);
    INFO(bformat("Alignment reads using %d threads starts") % param.threads);

    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 8))
    for (size_t j = 0; j < reads.size(); ++j) {
        CharString read_id = read_ids[j];

        // temporary solution - we crop three nucleotides at the end of string
        if (length(reads[j]) <= param.num_cropped_nucls) {
            // Discard so short read
            output_isok[j] = false;
            continue;
        }

        Dna5String read = prefix(reads[j], length(reads[j]) - param.num_cropped_nucls);
        Dna5String read_rc = read;
        reverseComplement(read_rc);

        auto result_pstrand = index.query(read, param.max_candidates);
        auto result_nstrand = index.query(read_rc, param.max_candidates);

        int pscore = (result_pstrand.size() > 0) ? result_pstrand[0].kp_coverage : 0;
        int nscore = (result_nstrand.size() > 0) ? result_nstrand[0].kp_coverage : 0;

        int strand = (pscore >= nscore) ? 1 : -1;
        auto result = (strand == 1) ? result_pstrand : result_nstrand;
        auto stranded_read = (strand == 1) ? read : read_rc;

        bool aligned = false;
        if (!result.empty()) { // If we found at least one alignment
            for (const auto &align : result) {
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

                    cout << " " << visualize_matches(align.path, length(v_reads[align.needle_index]), length(stranded_read)) << endl;
                    if (!jresult.empty()) {
                        cout << "\tJ-genes:" << endl;
                        for (const auto &jalign : jresult) {
                            cout << bformat("\tk+-coverage %3%; %1%:%2%; J-gene: %4%")
                                % (jalign.start+1 + end_of_v) % (jalign.finish + end_of_v)
                                % jalign.kp_coverage          % j_ids[jalign.needle_index];

                            cout << " " << visualize_matches(jalign.path, length(j_reads[jalign.needle_index]), length(suff)) << endl;
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

                    {
                        const auto &first_jalign = jalign.path[0];
                        const auto &last_jalign = jalign.path[jalign.path.size() - 1];

                        bformat bf("%s, %d, %d, %1.2f, %s, %d, %d, %1.2f, %s");
                        bf % read_id
                           % (align.start+1)             % end_of_v
                           % align.score                 % toCString(v_ids[align.needle_index])
                           % (first_jalign.read_pos + 1 + end_of_v) % (last_jalign.read_pos + last_jalign.length + end_of_v)
                           % jalign.score                % toCString(j_ids[jalign.needle_index]);

                        add_info_strings[j] = bf.str();

                        reads[j] = cropped_read;
                        output_isok[j] = true;
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
    add_info << bformat("%s, %s, %s, %s, %s, %s, %s, %s, %s\n")
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
