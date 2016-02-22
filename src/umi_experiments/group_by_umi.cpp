#include <ostream>
#include <iostream>
#include <logger/log_writers.hpp>
#include <seqan/seq_io.h>
#include <segfault_handler.hpp>
#include "../fast_ig_tools/banded_half_smith_waterman.hpp"
#include "umi_utils.hpp"
#include "utils.hpp"

bool parse_cmdline(int argc, char **argv, std::string& input_file, std::string& output_dir) {
    if (argc != 3) {
        std::cout << "Usage: <input file> <output dir>" << std::endl;
        return false;
    }
    input_file = argv[1];
    output_dir = argv[2];
    return true;
}

class ReadGroup {
public:
    ReadGroup(const seqan::Dna5String& first/*, seqan::DnaQString& qual*/) : center_(first)/*, min_length_(length(first)), nt_representation_(length(first))*/ {
        representatives_.push_back(first);
//        for (size_t i = 0; i < nt_representation_.size(); i ++) {
//            nt_representation_[i] = vector<char>(4);
//        }
    }

    ReadGroup() = default;

    bool TryAdd(const seqan::Dna5String& candidate/*, seqan::DnaQString& quality*/) {
        if (length(center_) == 0) {
            center_ = candidate;
            return true;
        }
        auto dist = [&](const seqan::Dna5String& s1, const seqan::Dna5String& s2) -> int {
            return -half_sw_banded(s1, s2, 0, -1, -1, [](int) -> int { return 0; }, ReadGroup::indels_);
        };
        if (dist(candidate, center_) > ReadGroup::radius_) {
            return false;
        }
        representatives_.push_back(candidate);
//        if (representatives_.size() > FIND_CENTER_EVERY_TIME_THRESHOLD && representatives_.size() % REEVALUATE_CENTER_INTERVAL != 0) {
//            return true;
//        }
        reevaluate_center(/*candidate, dist*/);
        return true;
    }

    size_t Size() const { return representatives_.size(); }

    std::vector<seqan::Dna5String> GetRepresentatives() const { return representatives_; }

    seqan::Dna5String GetCenter() const { return center_; }

    static void SetRadius(int radius) { radius_ = radius; }

    static void SetTau(int indels) { indels_ = indels; }

private:
    static const size_t FIND_CENTER_EVERY_TIME_THRESHOLD = 5;
    static const size_t REEVALUATE_CENTER_INTERVAL = 5;

    static int radius_;
    static int indels_;

    seqan::Dna5String center_;
    std::vector<seqan::Dna5String> representatives_;

//    size_t min_length_;
//    vector<vector<char>> nt_representation_;

    void reevaluate_center(/*const seqan::Dna5String& candidate, function<int(const seqan::Dna5String&, const seqan::Dna5String&)> dist*/) {
        auto len = length(center_);
        for (size_t i = 0; i < len; i ++) {
            std::vector<size_t> cnt(4);
            for (auto repr : representatives_) {
                if (length(repr) > i) {
                    cnt[repr[i]]++;
                }
            }
            for (size_t j = 0; j < 4; j ++) {
                if (cnt[j] > cnt[center_[i]]) {
                    center_[i] = j;
                }
            }
        }

//        min_length_ = min(min_length_, length(candidate));
//        for (size_t i = 0; i < min_length_; i ++) {
//            nt_representation_[i][candidate[i]] ++;
//            for (size_t nt = 0; nt < nt_representation_[i].size(); nt ++) {
//                if (nt_representation_[nt] > nt_representation_[center_[i]]) {
//                    center_[i] = nt;
//                }
//            }
//        }
    }
};

int ReadGroup::radius_ = 0;
int ReadGroup::indels_ = 0;

class UmiReadSet {
public:
    UmiReadSet() {}

    void AddRead(const seqan::Dna5String& read/*, seqan::DnaQString& qual*/) {
        for (size_t i = 0; i < readGroups_.size(); i ++) {
            if (readGroups_[i].TryAdd(read/*, qual*/)) {
//                read_to_group_[read] = i;
                return;
            }
        }

        readGroups_.emplace_back(read);
//        read_to_group_[read] = readGroups_.size() - 1;
    }

    const std::vector<ReadGroup>& GetReadGroups() const { return readGroups_; }

private:
    std::vector<ReadGroup> readGroups_;
//    unordered_map<seqan::Dna5String, size_t> read_to_group_;
};

void group_by_umi(std::vector<seqan::Dna5String>& input_umi, std::vector<seqan::DnaQString>& input_umi_qual,
                  std::vector<seqan::Dna5String>& input_reads, std::vector<seqan::DnaQString>& input_qual,
                  std::unordered_map<Umi, UmiReadSet>& umi_to_reads) {
    VERIFY_MSG(input_umi.size() == input_umi_qual.size() && input_umi.size() == input_reads.size() && input_umi.size() == input_qual.size(), "Bad info read");
    size_t mean_read_length = 0;
    for (auto& read: input_reads) {
        mean_read_length += length(read);
    }
    mean_read_length /= input_reads.size();

    for (int indels = 0; indels <= 10; indels ++) {
        ReadGroup::SetTau(indels);
        for (int radius = /*static_cast<int>(mean_read_length) / 4*/static_cast<int>(mean_read_length), ridx = 0; /*ridx < 1*/radius > 0; radius /= 2, ridx ++) {
            ReadGroup::SetRadius(radius);
            INFO("Calculating for indels = " << indels << " and radius " << radius);
            umi_to_reads.clear();
            for (size_t i = 0; i < input_umi.size(); i ++) {
                umi_to_reads[Umi(input_umi[i])].AddRead(input_reads[i]);
            }

            INFO("Counting");
            size_t max_group_size = 0;
            size_t total_group_size = 0;
            size_t group_count = 0;
            std::vector<size_t> group_sizes(10);
            size_t max_groups_in_umi = 0;
            size_t min_center_length = std::numeric_limits<size_t>::max();
            for (auto& entry : umi_to_reads) {
                auto read_groups = entry.second.GetReadGroups();
                size_t umi_size = 0;
                for (auto& group : read_groups) {
                    max_group_size = std::max(max_group_size, group.Size());
                    total_group_size += group.Size();
                    umi_size += group.Size();
                    group_count ++;
                    while (group.Size() >= group_sizes.size()) {
                        group_sizes.push_back(0);
                    }
                    group_sizes[group.Size()] ++;
                    min_center_length = std::min(min_center_length, length(group.GetCenter()));
//                    for (auto read : group.GetRepresentatives()){
//                        cout << read << endl;
//                    }
//                    cout << "------------" << endl;
//                    cout << group.GetCenter() << endl;
//                    if (group.Size() >= 2) return;
                }
//                cout << "-----------------------------" << endl;
                max_groups_in_umi = std::max(max_groups_in_umi, read_groups.size());
//                if (umi_size >= 2) return;
            }
            VERIFY_MSG(input_umi.size() == total_group_size, "Lost or acquired extra reads. Should be " << input_umi.size() << ", but got " << total_group_size);
            INFO("Tau " << indels << ", radius " << radius << ": unique UMIs " << umi_to_reads.size() << " max groups in UMI " << max_groups_in_umi <<
                 " max group size " << max_group_size << ", mean group size " << (static_cast<double>(total_group_size) / static_cast<double>(group_count)) <<
                 " min center length " << min_center_length);
            std::stringstream ss;
            for (size_t size = 1; size <= max_group_size; size ++) {
                ss << boost::format("%d,") % group_sizes[size];
            }
            INFO(ss.str());
        }
    }
}

int main(int argc, char** argv) {
    segfault_handler sh;
    create_console_logger();

    std::string input_file;
    std::string output_dir;
    if (!parse_cmdline(argc, argv, input_file, output_dir)) return 0;

    INFO("Reading fastq");
    seqan::SeqFileIn seqFileIn_input(input_file.c_str());
    std::vector<seqan::CharString> input_ids;
    std::vector<seqan::Dna5String> input_reads;
    std::vector<seqan::DnaQString> input_qual;
    readRecords(input_ids, input_reads, input_qual, seqFileIn_input);

    INFO("Extracting UMI data");
    std::vector<seqan::Dna5String> input_umi;
    std::vector<seqan::DnaQString> input_umi_qual;
    extract_barcodes_from_read_ids(input_ids, input_umi, input_umi_qual);

    INFO("Grouping by UMI");
    std::unordered_map<Umi, UmiReadSet> umi_to_reads;
    group_by_umi(input_umi, input_umi_qual, input_reads, input_qual, umi_to_reads);

    // construct umi graph split using read list data
}
