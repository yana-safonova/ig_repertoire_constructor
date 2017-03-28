#include <cassert>
#include <algorithm>

#include <unordered_map>
#include <build_info.hpp>

#include <iostream>
#include <sstream>
using std::cout;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/algorithm/string.hpp>

#include "fast_ig_tools.hpp"
#include "ig_final_alignment.hpp"
#include "ig_matcher.hpp"
#include "utils.hpp"

#include <seqan/seq_io.h>
#undef NDEBUG
#include <cassert>

using seqan::Dna5String;
using seqan::SeqFileIn;
using seqan::SeqFileOut;
using seqan::CharString;


std::unordered_map<std::string, std::string> read_rcm_file(const std::string &file_name) {
    std::ifstream rcm(file_name.c_str());

    std::unordered_map<std::string, std::string> result;

    std::string id, target;
    std::string line;
    while (std::getline(rcm, line)) {
        std::vector<std::string> strs;
        boost::split(strs, line, boost::is_any_of("\t"));
        for (auto &s : strs) {
            boost::trim(s);
        }
        if (strs.size() == 2 && strs[0] != "" && strs[1] != "") {
            result[strs[0]] = strs[1];
        }
    }

    return result;
}


template<typename T = seqan::Dna5>
void split_component(const std::vector<seqan::String<T>> &reads,
                     const std::vector<size_t> &indices,
                     std::vector<std::pair<seqan::String<T>, std::vector<size_t>>> &out,
                     size_t max_votes = 1,
                     bool discard = false,
                     bool recursive = true,
                     bool flu = true) {
    if (!max_votes) {
        max_votes = std::numeric_limits<size_t>::max() / 2;
    }

    if (indices.size() == 0) {
        return;
    }

    if (indices.size() == 1) {
        out.push_back({ reads[indices[0]], indices });
        return;
    }

    using namespace seqan;

    String<ProfileChar<T>> profile;

    size_t len = 0;
    for (size_t i : indices) {
        len = std::max(len, length(reads[i]));
    }

    resize(profile, len);

    for (size_t i : indices) {
        const auto &read = reads[i];
        for (size_t j = 0; j < length(read); ++j) {
            profile[j].count[ordValue(read[j])] += 1;
        }
    }

    // Find secondary votes
    struct PositionVote {
        size_t majory_votes;
        size_t majory_letter;
        size_t secondary_votes;
        size_t secondary_letter;
        size_t position;
        bool operator<(const PositionVote &b) const {
            return secondary_votes < b.secondary_votes;
        }
    };

    // INFO("Splitting component size=" << indices.size() << " len=" << len);
    size_t min_len = length(reads[indices[0]]);
    for (size_t i : indices) {
        min_len = std::min(min_len, length(reads[i]));
    }

    std::vector<PositionVote> secondary_votes;
    for (size_t j = 0; j < min_len; ++j) {
        std::vector<std::pair<size_t, size_t>> v;
        for (size_t k = 0; k < 4; ++k) {
            v.push_back({ profile[j].count[k], k });
        }

        // Use nth element here???
        std::sort(v.rbegin(), v.rend());
        secondary_votes.push_back({ v[0].first, v[0].second, v[1].first, v[1].second, j });
    }

    auto maximal_mismatch = *std::max_element(secondary_votes.cbegin(), secondary_votes.cend());
    VERIFY(maximal_mismatch.majory_votes >= maximal_mismatch.secondary_votes);

    TRACE("VOTES: " << maximal_mismatch.majory_votes << "/" << maximal_mismatch.secondary_votes << " POSITION: " << maximal_mismatch.position);
    bool do_split = false;
    auto mmsv = maximal_mismatch.secondary_votes;

    if (flu) {
        do_split = -0.0064174097073423269 * static_cast<double>(indices.size()) + 0.79633984048973583 * static_cast<double>(mmsv) - 4.3364230321953841 > 0;
    } else {
        do_split = mmsv >= max_votes;
    }
    if (indices.size() <= 5) {
        do_split = false;
    }

    if (max_votes > indices.size()) {
        do_split = false;
    }

    if (! do_split) {
        seqan::String<T> consensus;
        for (size_t i = 0; i < length(profile); ++i) {
            size_t idx = getMaxIndex(profile[i]);
            if (idx < ValueSize<T>::VALUE) {  // is not gap  TODO Check it!!
                appendValue(consensus, T(getMaxIndex(profile[i])));
            }
        }

        out.push_back({ consensus, indices });
        return;
    }

    std::vector<size_t> indices_majory, indices_secondary, indices_other;
    for (size_t i : indices) {
        if (seqan::ordValue(reads[i][maximal_mismatch.position]) == maximal_mismatch.majory_letter) {
            indices_majory.push_back(i);
        } else if (seqan::ordValue(reads[i][maximal_mismatch.position]) == maximal_mismatch.secondary_letter) {
            indices_secondary.push_back(i);
        } else {
            indices_other.push_back(i);
        }
    }

    VERIFY(indices_majory.size() == maximal_mismatch.majory_votes);
    VERIFY(indices_secondary.size() == maximal_mismatch.secondary_votes);

    auto majory_consensus = consensus_hamming(reads, indices_majory);
    auto secondary_consensus = consensus_hamming(reads, indices_secondary);

    for (size_t i : indices_other) {
        auto dist_majory = hamming_rtrim(reads[i], majory_consensus);
        auto dist_secondary = hamming_rtrim(reads[i], secondary_consensus);

        if (dist_majory <= dist_secondary) {
            indices_majory.push_back(i);
        } else {
            indices_secondary.push_back(i);
        }
    }

    VERIFY(indices_majory.size() + indices_secondary.size() == indices.size());
    VERIFY(indices_majory.size() <= indices.size());

    INFO("Component splitted " << indices_majory.size() << " + " << indices_secondary.size());

    if (!recursive) {
        max_votes = 0;
    }

    split_component(reads, indices_majory, out, max_votes, discard, flu);

    if (discard) {
        for (size_t index : indices_secondary) {
            out.push_back({ reads[index], { index } });
        }
    } else {
        VERIFY(indices_secondary.size() < indices.size());
        split_component(reads, indices_secondary, out, max_votes, discard, flu);
    }
}


template<typename T = seqan::Dna5>
std::vector<std::pair<seqan::String<T>, std::vector<size_t>>> split_component(const std::vector<seqan::String<T>> &reads,
                                                                              const std::vector<size_t> &indices,
                                                                              size_t max_votes = 0,
                                                                              bool discard = false,
                                                                              bool recursive = true,
                                                                              bool flu = true) {
    if (!max_votes) {
        max_votes = std::numeric_limits<size_t>::max() / 2;
    }

    std::vector<std::pair<seqan::String<T>, std::vector<size_t>>> result;
    split_component(reads, indices, result, max_votes, discard, recursive, flu);

    return result;
}


int main(int argc, char **argv) {
    segfault_handler sh;
    perf_counter pc;
    create_console_logger("");

    INFO("Command line: " << join_cmd_line(argc, argv));

    int nthreads = 4;
    std::string reads_file;
    std::string output_file;
    std::string rcm_file;
    std::string output_rcm_file;
    std::string config_file;
    size_t max_votes = std::numeric_limits<size_t>::max() / 2;
    bool discard = false;
    bool recursive = true;
    bool flu = false;
    bool allow_unassigned = false;

    // Parse cmd-line arguments
    try {
        // Declare a group of options that will be
        // allowed only on command line
        po::options_description generic("Generic options");
        generic.add_options()
            ("version,v", "print version string")
            ("help,h", "produce help message")
            ("config,c", po::value<std::string>(&config_file)->default_value(config_file),
             "name of a file of a configuration")
            ("input-file,i", po::value<std::string>(&reads_file),
             "name of the input file (FASTA|FASTQ)")
            ("output-file,o", po::value<std::string>(&output_file)->default_value(output_file),
             "output file for final repertoire")
            ("rcm-file,R", po::value<std::string>(&rcm_file)->default_value(rcm_file),
             "input RCM-file")
            ("output-rcm-file,M", po::value<std::string>(&output_rcm_file)->default_value(output_rcm_file),
             "output RCM-file");

        // Declare a group of options that will be
        // allowed both on command line and in
        // config file
        po::options_description config("Configuration");
        config.add_options()
            ("threads,t", po::value<int>(&nthreads)->default_value(nthreads),
             "the number of parallel threads")
            ("max-votes,V", po::value<size_t>(&max_votes)->default_value(max_votes),
             "max secondary votes threshold")
            ("discard,D", po::value<bool>(&discard)->default_value(discard),
             "whether to discard secondary votes")
            ("recursive,C", po::value<bool>(&recursive)->default_value(recursive),
             "whether to perform recursive splitting")
            ("flu,F", po::value<bool>(&flu)->default_value(flu),
             "Use FLU preset")
            ("allow-unassigned,A", po::value<bool>(&allow_unassigned)->default_value(allow_unassigned),
             "Allow unassigned reads in RCM")
            ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        po::options_description hidden("Hidden options");
        hidden.add_options()
            ("help-hidden", "show all options, including developers' ones")
            ;

        po::options_description cmdline_options("All command line options");
        cmdline_options.add(generic).add(config).add(hidden);

        po::options_description config_file_options;
        config_file_options.add(config).add(hidden);

        po::options_description visible("Allowed options");
        visible.add(generic).add(config);

        po::positional_options_description p;
        p.add("input-file", 1);

        po::variables_map vm;
        store(po::command_line_parser(argc, argv).
              options(cmdline_options).run(), vm);
        notify(vm);

        if (vm.count("help-hidden")) {
            cout << cmdline_options << std::endl;
            return 0;
        }

        if (vm.count("help")) {
            cout << visible << std::endl;
            return 0;
        }

        if (vm.count("version")) {
            cout << bformat("IG Component Splitter, part of IgReC version %s; git version: %s") % build_info::version % build_info::git_hash7 << std::endl;
            return 0;
        }

        if (vm.count("config-file")) {
            std::string config_file = vm["config-file"].as<std::string>();

            std::ifstream ifs(config_file.c_str());
            if (!ifs) {
                cout << "can not open config file: " << config_file << "\n";
                return 0;
            } else {
                store(parse_config_file(ifs, config_file_options), vm);
                // reparse cmd line again for update config defaults
                store(po::command_line_parser(argc, argv).
                      options(cmdline_options).positional(p).run(), vm);
            }
        }

        try {
            notify(vm);
        } catch (po::error &e) {
            cout << "Parser error: " << e.what() << std::endl;
        }
    } catch(std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    INFO("Command line: " << join_cmd_line(argc, argv));
    INFO("Input files: " << reads_file << ", " << rcm_file);

    std::vector<Dna5String> input_reads;
    std::vector<CharString> input_ids;
    SeqFileIn seqFileIn_input(reads_file.c_str());

    INFO("Reading input reads starts");
    readRecords(input_ids, input_reads, seqFileIn_input);
    INFO(input_reads.size() << " reads were extracted from " << reads_file);

    std::vector<size_t> component_indices;
    component_indices.resize(input_reads.size());

    INFO("Reading read-cluster map starts");
    const auto rcm = read_rcm_file(rcm_file);

    std::unordered_map<std::string, std::vector<size_t>> comp2readnum;

    size_t assigned_reads = 0;
    for (size_t i = 0; i < input_reads.size(); ++i) {
        std::string id = toCString(input_ids[i]);
        auto pcomponent = rcm.find(id);
        if (pcomponent != rcm.end()) {
            comp2readnum[pcomponent->second].push_back(i);
            ++assigned_reads;
        } else {
            TRACE("Read " << id <<  " not found in RCM " << rcm_file);
        }
    }

    INFO(assigned_reads << " over " << input_reads.size() << " reads assigned");
    if (assigned_reads < input_reads.size() && !allow_unassigned) {
        ERROR(input_reads.size() - assigned_reads << " unassigned reads in RCM " << rcm_file);
        ERROR("Unassigned reads are not allowed");
        ERROR("Pass option '--allow-unassigned=true' to allow");
        return 1;
    }

    INFO(comp2readnum.size() << " clusters were extracted from " << rcm_file);

    size_t max_component_size = 0;
    for (const auto &kv : comp2readnum) {
        max_component_size = std::max(max_component_size, kv.second.size());
    }
    INFO(bformat("Size of maximal cluster: %d") % max_component_size);

    omp_set_num_threads(nthreads);
    INFO(bformat("Computation of consensus using %d threads starts") % nthreads);
    INFO("Saving results");

    SeqFileOut seqFileOut_output(output_file.c_str());

    std::ofstream out_rcm(output_rcm_file.c_str());


    std::vector<std::pair<std::string, std::vector<size_t>>> comp2readnum_sorted(comp2readnum.cbegin(), comp2readnum.cend());
    std::sort(comp2readnum_sorted.begin(), comp2readnum_sorted.end());

    // TODO Use multithreading
    for (const auto &kv : comp2readnum_sorted) {
        const auto &comp = kv.first;
        const auto &indices = kv.second;
        auto result = split_component(input_reads, indices, max_votes, discard, recursive, flu);
        for (size_t i = 0; i < result.size(); ++i) {
            std::stringstream ss(comp);
            if (result.size() > 1) {
                ss << "X" << i;
            }
            std::string cluster_id = ss.str();

            bformat fmt("cluster___%s___size___%d");
            fmt % cluster_id % result[i].second.size();
            std::string id = fmt.str();

            seqan::writeRecord(seqFileOut_output, id, result[i].first);
            for (size_t read_index : result[i].second) {
                std::string read_id = seqan::toCString(input_ids[read_index]);
                out_rcm << read_id << "\t" << cluster_id << "\n";
            }
        }
    }

    INFO("Final repertoire was written to " << output_file);
    INFO("Final RCM was written to " << output_rcm_file);

    INFO("Running time: " << running_time_format(pc));

    return 0;
}

// vim: ts=4:sw=4
