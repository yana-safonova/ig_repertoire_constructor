#define SEQAN_HAS_ZLIB 1 // TODO Set it using cmake
#define SEQAN_HAS_BZLIB 1

#include <seqan/find.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/modifier.h>
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


using namespace seqan;
using std::vector;
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::map;
using std::make_pair;
using std::make_tuple;
using bformat = boost::format;

#include <boost/program_options.hpp>
namespace po = boost::program_options;


// input vector(pair(needle_pos, read_pos))
// output vector(tuple(needle_pos, read_pos, length))
struct Match {
  int needle_pos;
  int read_pos;
  unsigned length;
};


int path_coverage_length(const std::vector<Match> &path) {
  int result = 0;

  for (const auto &match : path) {
    result += match.length;
  }

  return result;
}


std::vector<Match> combine_sequential_kmer_matches(std::vector<std::pair<int, int>> matches,
                                                   unsigned K) {
  // Convert to (shift, read_pos)
  for (auto &match : matches) {
    match.first = match.second - match.first;
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
  using PatternIndex = Pattern<StringSet<Dna5String>, AhoCorasick>;
public:
  KmerIndex(StringSet<Dna5String> queries, int K,
            int max_global_gap, int left_uncoverage_limit, int right_uncoverage_limit,
            int max_local_insertions, int max_local_deletions, int min_k_coverage) : max_local_insertions{max_local_insertions},
                                                                                     max_local_deletions{max_local_deletions},
                                                                                     min_k_coverage{min_k_coverage} {
    this->queries = queries; // FIXME
    this->K = K;
    this->left_uncoverage_limit = left_uncoverage_limit;
    this->right_uncoverage_limit = right_uncoverage_limit;
    this->max_global_gap = max_global_gap;

    for (size_t j = 0; j < length(queries); ++j) {
      for (size_t start = 0; start + K <= length(queries[j]); start += 1) {
        kmer2needle[infix(queries[j], start, start + K)].push_back(std::make_pair(j, start));
      }
    }

    // Reverse for LIS search
    // Very important step!
    // We don't use LIS more...
    // for (auto &p : kmer2needle) {
    //   std::reverse(p.second.begin(), p.second.end());
    // }

    most_pop_kmer_uses = 0;
    for (const auto &e : kmer2needle) {
      most_pop_kmer_uses = std::max<int>(most_pop_kmer_uses, e.second.size());
    }

    StringSet<Dna5String> kmers;
    for (const auto &e : kmer2needle) {
      appendValue(kmers, e.first);
    }

    // Define the Aho-Corasick pattern over the queries with the preprocessing
    // data structure.
    // Better to use make_unique
    pattern = std::move(std::unique_ptr<PatternIndex>(new PatternIndex(kmers)));
    pmtx = std::move(std::unique_ptr<std::mutex>(new std::mutex));
  }

  std::unordered_map<size_t, std::vector<std::pair<int, int>> > Needle2matches(Dna5String read) {
    std::unordered_map<size_t, std::vector<std::pair<int, int>> > needle2matches;

    Finder<Dna5String> finder(read);  // new finder for each seq
    // Lock it
    std::lock_guard<std::mutex> lck(*this->pmtx);

    while (find(finder, *pattern)) {
      auto kmer = infix(finder);

      for (const auto &p : kmer2needle[kmer]) {
        size_t needle_index = p.first;
        int kmer_pos_in_read = position(finder);
        int kmer_pos_in_needle = p.second;
        int shift = kmer_pos_in_read - kmer_pos_in_needle;

        // We make these limits less strict because of possibility of indels
        int shift_min = -left_uncoverage_limit - max_global_gap;
        int shift_max = int(length(read)) - int(length(queries[needle_index])) + right_uncoverage_limit + max_global_gap;

        if (shift >= shift_min && shift <= shift_max) {
          needle2matches[needle_index].push_back(make_pair(kmer_pos_in_needle, kmer_pos_in_read));
        }
      }
    }

    return needle2matches;
  }

  struct Alignment {
    int kp_coverage;
    std::vector<Match> path;
    int start, finish;
    size_t needle_index;
    int overlap_length;

    bool operator< (const Alignment& b) const {
      return this->kp_coverage < b.kp_coverage;
    }
  };

  std::vector<Alignment> query_unordered(Dna5String read) {
    std::vector<Alignment> liss;
    auto needle2matches = this->Needle2matches(read);

    if (!needle2matches.empty()) { // if we found at least one proper k-mer match
      for (const auto &p : needle2matches) {
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
        align.path = std::move(path);
        align.start = start;
        align.finish = finish;
        align.needle_index = needle_index;
        align.overlap_length = needle_overlap_length;

        liss.push_back(std::move(align));
      }
    }

    return liss;
  }

  std::vector<Alignment> query(Dna5String read, int limit) {
    auto liss = query_unordered(read);

    using ctuple_type = decltype(*liss.cbegin());

    auto score_function = [](ctuple_type a) { return a.kp_coverage; };
    auto comp = [&](ctuple_type a, ctuple_type b) -> bool { return score_function(a) > score_function(b); };

    // Return top <limit>
    std::nth_element(liss.begin(), liss.begin() + limit, liss.end(), comp);
    liss.resize(std::min<int>(liss.size(), limit));
    std::sort(liss.begin(), liss.end(), comp);

    return liss;
  }

private:
  StringSet<Dna5String> queries;
  int max_global_gap;
  int left_uncoverage_limit, right_uncoverage_limit;
  int max_local_insertions;
  int max_local_deletions;
  int min_k_coverage;
  int most_pop_kmer_uses = 0;
  int K;
  std::unique_ptr<PatternIndex> pattern;
  map<Dna5String, vector<std::pair<int, int> > > kmer2needle;
  std::unique_ptr<std::mutex> pmtx;
};


int main(int argc, char **argv) {
  int K = 7; // anchor length
  int word_size_j = 5;
  int left_uncoverage_limit = 16;
  int right_uncoverage_limit = 0; // We have to cover (D)J region too. Maybe it should be even negative
  std::string input_file = "", organism = "human";
  int max_local_deletions = 9;
  int max_local_insertions = 9;
  int max_global_gap = 18;
  int min_k_coverage = 50;
  int min_k_coverage_j = 9;
  int max_candidates = 3;
  int max_candidates_j = 3;
  std::string chain = "heavy";
  std::string db_directory = "./germline";
  std::string output_filename = "cropped.fa";
  std::string bad_output_filename = "bad.fa";
  std::string add_info_filename = "add_info.csv";
  int threads = 4;
  bool silent = true;
  bool fill_prefix_by_germline = true;

  // Parse cmd-line arguments
  try {
    std::string config_file;

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
      ("silent,S", "suppose output for each query (default)")
      ("no-silent,V", "produce info output for each query")
      ("chain,C", po::value<std::string>(&chain)->default_value(chain),
       "IG chain ('heavy', 'lambda' or 'kappa')")
      ("db-directory", po::value<std::string>(&db_directory)->default_value(db_directory),
       "directory with germline database")
      ("output,o", po::value<std::string>(&output_filename)->default_value(output_filename),
       "output file for <<good>> reads ([right-stranded], filtered and cropped)")
      ("add-info,a", po::value<std::string>(&add_info_filename)->default_value(add_info_filename),
       "output file for additional aligment info (in CSV format)")
      ("bad-output,b", po::value<std::string>(&bad_output_filename)->default_value(bad_output_filename),
       "output file for <<bad>> reads")
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
        cout << "can not open config file: " << config_file << "\n";
        return 0;
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
      return 0;
    }

    if (vm.count("help") || !vm.count("input-file")) { // TODO Process required arguments by the proper way
      cout << visible << "\n";
      return 0;
    }

    if (vm.count("version")) {
      cout << "<Some cool name> version 0.1" << vm.count("version") << std::endl;
      return 0;
    }

    if (vm.count("input-file")) {
      cout << "Input file is: "
        << vm["input-file"].as<std::string>() << "\n";
    }

    if (vm.count("silent")) {
      silent = true;
    } else if (vm.count("no-silent")) {
      silent = false;
    }

    if (vm.count("no-fill-prefix-by-germline")) {
      fill_prefix_by_germline = false;
    }

    if (vm.count("fill-prefix-by-germline")) {
      fill_prefix_by_germline = true;
    }

    cout << "K = " << K << endl;
    cout << bformat("Output files are: %s %s %s") % output_filename % bad_output_filename % add_info_filename << std::endl;
  } catch(std::exception& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  vector<CharString> v_ids;
  StringSet<Dna5String> v_reads;

  std::string chain_letter;
  if (chain == "heavy") {
    chain_letter = "H";
  } else if (chain == "lambda") {
    chain_letter = "L";
  } else if (chain == "kappa") {
    chain_letter = "K";
  } else {
    std::cerr << "chain should be 'heavy', 'lambda' or 'kappa'" << std::endl;
    return 1;
  }

  std::string vgenes_filename = db_directory + "/" + organism + "/IG" + chain_letter + "V.fa";
  std::string jgenes_filename = db_directory + "/" + organism + "/IG" + chain_letter + "J.fa";
  SeqFileIn seqFileIn_v_genes(vgenes_filename.c_str());
  SeqFileIn seqFileIn_j_genes(jgenes_filename.c_str());

  readRecords(v_ids, v_reads, seqFileIn_v_genes);
  vector<CharString> j_ids;
  StringSet<Dna5String> j_reads;
  readRecords(j_ids, j_reads, seqFileIn_j_genes);

  seqan::SeqFileOut cropped_reads_seqFile(output_filename.c_str());
  seqan::SeqFileOut bad_reads_seqFile(bad_output_filename.c_str());
  std::ofstream add_info(add_info_filename.c_str());
  add_info << bformat("%s, %s, %s, %s, %s, %s, %s, %s\n")
    % "Vstart" % "Vend"
    % "Vscore" % "Vindex"
    % "Jstart" % "Jend"
    % "Jscore" % "Jindex";

  seqan::SeqFileIn seqFileIn_reads(input_file.c_str());

  std::mutex read_mtx, write_mtx, stdout_mtx;

  auto work = [&](KmerIndex index, KmerIndex j_index) -> void {
    while (true) {
    CharString read_id;
    Dna5String read;

    {
      std::lock_guard<std::mutex> lck(read_mtx);
      if (!atEnd(seqFileIn_reads)) {
        readRecord(read_id, read, seqFileIn_reads);
      } else {
        return;
      }
    }

    Dna5String read_rc = read;
    reverseComplement(read_rc);

    auto result_pstrand = index.query(read, max_candidates);
    auto result_nstrand = index.query(read_rc, max_candidates);

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
        auto jresult = j_index.query(suff, max_candidates_j);

        if (!silent) {
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
            if (fill_prefix_by_germline) {
              cropped_read = prefix(v_reads[align.needle_index], std::abs(align.start)); // Fill by germline
            } else {
              cropped_read = std::string(std::abs(align.start), 'N');
            }
            cropped_read += stranded_read;
          }
          // if (align.finish > length(stranded_read)) {
          //   cropped_read += std::string(align.finish - length(stranded_read), 'N');
          // }

          const auto &jalign = *jresult.cbegin();

          {
            const auto &first_jalign = jalign.path[0];
            const auto &last_jalign = jalign.path[jalign.path.size() - 1];

            std::lock_guard<std::mutex> lck(write_mtx);
            add_info << bformat("%d, %d, %d, %d, %d, %d, %d, %d\n")
              % (align.start+1)             % end_of_v
              % align.kp_coverage           % align.needle_index
              % (first_jalign.read_pos + 1 + end_of_v) % (last_jalign.read_pos + last_jalign.length + end_of_v)
              % jalign.kp_coverage          % jalign.needle_index;

            seqan::writeRecord(cropped_reads_seqFile, read_id, cropped_read);
          }

          aligned = true;
        }
      }
    } else if (!silent) {
      std::lock_guard<std::mutex> lck(stdout_mtx);

      cout << "Query: " << read_id << endl;
      cout << bformat("Strand: %s (best V-gene k+-score: +%d/-%d)")
        % ((strand == 1) ? "+" : "-") % pscore % nscore <<  endl; // '?' has less priority
      cout << "V-genes not found" << endl;
      cout << endl;
    }

    if (!aligned) { // Save as bad read
      std::lock_guard<std::mutex> lck(write_mtx);
      seqan::writeRecord(bad_reads_seqFile, read_id, read); // Write original read, not stranded
    }
  }
  };

  std::vector<std::thread> workers;
  for (int i = 0; i < threads; ++i) {
    workers.push_back(std::thread(work,
                                  KmerIndex(v_reads, K,
                                            max_global_gap, left_uncoverage_limit, right_uncoverage_limit,
                                            max_local_insertions, max_local_deletions, min_k_coverage),
                                  KmerIndex(j_reads, word_size_j,
                                            max_global_gap, 100000, 10000,
                                            max_local_insertions, max_local_deletions, min_k_coverage_j)));
  }

  for (auto &t : workers) {
    t.join();
  }

  return 0;
}
