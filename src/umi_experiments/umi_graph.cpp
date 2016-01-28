#include <bits/stringfwd.h>
#include <iostream>
#include <logger/log_writers.hpp>

#include <seqan/seq_io.h>
#include <segfault_handler.hpp>
#include "../ig_tools/utils/string_tools.hpp"
#include "../fast_ig_tools/banded_half_smith_waterman.hpp"

using seqan::Dna5String;
using seqan::DnaQString;
using seqan::SeqFileIn;
using seqan::CharString;

void create_console_logger() {
    using namespace logging;
    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

bool parse_cmdline(int argc, char **argv, std::string& file) {
    if (argc != 2) {
        std::cout << "Usage: <input file>" << std::endl;
        return false;
    }
    file = argv[1];
    return true;
}

void extract_umi(const std::vector<CharString>& ids, std::vector<Dna5String>& umis, std::vector<DnaQString>& umi_quals) {
    for (auto& id : ids) {
        auto id_string = seqan::toCString(id);
        auto split_by_umi = split(id_string, "UMI");
        VERIFY_MSG(split_by_umi.size() == 2, "Either no UMI info, or too much of it in meta: " << id_string);
        auto umi_info = split_by_umi[1].substr(1);
        auto colon = umi_info.find(":");
        VERIFY_MSG(colon != string::npos, "Can't parse UMI info: " << umi_info);
        auto umi = umi_info.substr(0, colon);
        auto qual = umi_info.substr(colon + 1);
        VERIFY_MSG(umi.length() == qual.length(), "UMI and its quality are of different lengths: " << umi_info);
        umis.push_back(umi);
        umi_quals.push_back(qual);
    }
}

class Umi {
public:
    Umi(const Dna5String& umi) : umi_(umi) {}

    bool operator==(const Umi &other) const { return umi_ == other.umi_; }

    Dna5String GetString() const { return umi_; }

private:
    const Dna5String umi_;
};

namespace std {
    template<>
    struct hash<Umi> {
        size_t operator()(const Umi& umi) const {
            size_t h = 0;
            for (auto c : umi.GetString()) {
                h = h * 31 + seqan::ordValue(c);
            }
            return h;
        }
    };
}

class ReadGroup {
public:
    ReadGroup(Dna5String& first/*, DnaQString& qual*/) : center_(first) {
        representatives_.push_back(first);
    }

    bool WouldLikeIn(Dna5String& candidate/*, DnaQString& quality*/, bool include) {
        auto dist = [&](const Dna5String& s1, const Dna5String& s2) -> int {
            return -half_sw_banded(s1, s2, 0, -1, -1, [](int) -> int { return 0; }, ReadGroup::tau_);
        };
        if (dist(candidate, center_) > ReadGroup::radius_) {
            return false;
        }
        if (include) {
            representatives_.push_back(candidate);
            Dna5String& best = center_;
            int best_dist = static_cast<int>(length(best)) * 2;
            for (auto& read : representatives_) {
                int max_dist = 0;
                for (auto& another : representatives_) {
                    max_dist = max(max_dist, dist(read, another));
                }
                if (max_dist < best_dist) {
                    best = read;
                    best_dist = max_dist;
                }
            }
            center_ = best;
        }
        return true;
    }

    size_t Size() const { return representatives_.size(); }

    static void SetRadius(int radius) { radius_ = radius; }

    static void SetTau(int tau) { tau_ = tau; }

private:
    static int radius_;
    static int tau_;

    Dna5String& center_;
    vector<Dna5String> representatives_;
};

int ReadGroup::radius_ = 0;
int ReadGroup::tau_ = 0;

class UmiReadSet {
public:
    UmiReadSet() {}

    void AddRead(Dna5String& read/*, DnaQString& qual*/) {
        for (auto& readGroup : readGroups_) {
            if (readGroup.WouldLikeIn(read/*, qual*/, true)) {
                return;
            }
        }
        readGroups_.push_back(ReadGroup(read/*, qual*/));
    }

    const vector<ReadGroup>& GetReadGroups() const { return readGroups_; }

private:
    vector<ReadGroup> readGroups_;
};

void group_by_umi(std::vector<Dna5String>& input_umi, std::vector<DnaQString>& input_umi_qual,
                  std::vector<Dna5String>& input_reads, std::vector<DnaQString>& input_qual,
                  unordered_map<Umi, UmiReadSet>& umi_to_reads) {
    VERIFY_MSG(input_umi.size() == input_umi_qual.size() && input_umi.size() == input_reads.size() && input_umi.size() == input_qual.size(), "Bad info read");
    size_t mean_read_length = 0;
    for (auto& read: input_reads) {
        mean_read_length += length(read);
    }
    mean_read_length /= input_reads.size();

    for (int tau = 0; tau <= 10; tau ++) {
        ReadGroup::SetTau(tau);
        for (int radius = static_cast<int>(mean_read_length), ridx = 0; ridx < 7; radius /= 2, ridx ++) {
            ReadGroup::SetRadius(radius);
            INFO("Calculating for tau = " << tau << " and radius " << radius);
            umi_to_reads.clear();
            for (size_t i = 0; i < input_umi.size(); i ++) {
//                size_t total_group_size = 0;
//                for (auto& entry : umi_to_reads) {
//                    for (auto& group : entry.second.GetReadGroups()) {
//                        total_group_size += group.Size();
//                    }
//                }
//                VERIFY_MSG(i == total_group_size, "Lost or acquired extra reads. Should be " << i << ", but got " << total_group_size);
//                if (i % 1000 == 0) {
//                    INFO("Adding " << i << " out of " << input_umi.size());
//                }
                umi_to_reads[input_umi[i]].AddRead(input_reads[i]/*, input_qual[i]*/);
            }

            INFO("Counting");
            size_t max_group_size = 0;
            size_t total_group_size = 0;
            size_t group_count = 0;
            for (auto& entry : umi_to_reads) {
                for (auto& group : entry.second.GetReadGroups()) {
                    max_group_size = max(max_group_size, group.Size());
                    total_group_size += group.Size();
                    group_count ++;
                }
            }
            VERIFY_MSG(input_umi.size() == total_group_size, "Lost or acquired extra reads. Should be " << input_umi.size() << ", but got " << total_group_size);
            INFO("Max group size " << max_group_size << ", mean group size " << (static_cast<double>(total_group_size) / static_cast<double>(group_count)));
        }
    }
}

int main(int argc, char** argv) {
    segfault_handler sh;
    create_console_logger();

    std::string input_file;
    if (!parse_cmdline(argc, argv, input_file)) return 0;

    INFO("Reading fastq");
    SeqFileIn seqFileIn_input(input_file.c_str());
    std::vector<CharString> input_ids;
    std::vector<Dna5String> input_reads;
    std::vector<DnaQString> input_qual;
    readRecords(input_ids, input_reads, input_qual, seqFileIn_input);

    INFO("Extracting UMI data");
    std::vector<Dna5String> input_umi;
    std::vector<DnaQString> input_umi_qual;
    extract_umi(input_ids, input_umi, input_umi_qual);

    INFO("Grouping by UMI");
    unordered_map<Umi, UmiReadSet> umi_to_reads;
    group_by_umi(input_umi, input_umi_qual, input_reads, input_qual, umi_to_reads);

    // construct umi graph split using read list data
}
