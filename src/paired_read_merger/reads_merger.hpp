#pragma once

#include <utils/fastq_reader.hpp>
#include <utils/sequence_tools.hpp>
#include <utils/string_tools.hpp>

#include "logger/log_writers.hpp"

struct merger_setting {
    size_t min_overlap;
    double max_mismatch_rate;
    bool simulated_mode;

    merger_setting() :
        min_overlap(50),
        max_mismatch_rate(.1),
        simulated_mode(false) { }

    void print() {
        INFO("Min overlap size: " << min_overlap);
        INFO("Max mismatch rate: " << max_mismatch_rate);
        if(simulated_mode)
            INFO("Simulated mode is ON");
    }
};

class SequenceMerger {
    merger_setting setting_;

    std::string reverse_quality_string(std::string qual) {
        for(size_t i = 0; i < qual.size() / 2; i++) {
            char tmp = qual[i];
            qual[i] = qual[qual.size() - i - 1];
            qual[qual.size() - i - 1] = tmp;
        }
        return qual;
    }

    std::pair<size_t, size_t> FindBestOverlap(std::string seq1, std::string seq2) {
        std::pair<size_t, size_t> best_overlap(size_t(-1), size_t(-1));
        for(size_t i = 0; i < seq1.size(); i++) {
            size_t overlap_size = std::min<size_t>(seq2.size(), seq1.size() - i);
            if(overlap_size >= setting_.min_overlap) {
                std::string subseq1 = seq1.substr(i, overlap_size);
                std::string subseq2 = seq2.substr(0, overlap_size);
                size_t dist = HammingDistance(subseq1, subseq2);
                double mism_rate = static_cast<double>(dist) / static_cast<double>(overlap_size);
                if(mism_rate <= setting_.max_mismatch_rate) {
                    if(best_overlap.first > dist) {
                        best_overlap.first = dist;
                        best_overlap.second = i;
                    }
                }
            }
        }
        return best_overlap;
    }

    std::pair<std::string, std::string> MergeQualifiedSeq(PairedFastqRead &paired_read) {
        std::string rc_right = reverse_complementary(paired_read.right_read.seq);
        std::string left = paired_read.left_read.seq;
        std::pair<size_t, size_t> overlap = FindBestOverlap(left, rc_right);
        if(overlap == std::make_pair(size_t(-1), size_t(-1)))
            return make_pair(std::string(), std::string());

        size_t overlap_size = std::min<size_t>(left.size() - overlap.second,
                rc_right.size());
        std::string overlap_seq = left.substr(overlap.second, overlap_size);
        std::string qual1 = paired_read.left_read.quality;
        std::string qual2 = reverse_quality_string(paired_read.right_read.quality);
        std::string overlap_qual = qual1.substr(overlap.second,
                overlap_size);

        if(!setting_.simulated_mode)
            for(size_t i = 0; i < overlap_size; i++)
                if(qual1[overlap.second + i] < qual2[i]) {
                    overlap_seq[i] = rc_right[i];
                    overlap_qual[i] = qual2[i];
                }

        std::string merged_seq = left.substr(0, overlap.second) + overlap_seq;
        std::string merged_qual = qual1.substr(0, overlap.second) + overlap_qual;
        if(overlap_size > rc_right.size()) {
            merged_seq = merged_seq + left.substr(overlap.second + overlap_size,
                    left.size() - overlap.second - overlap_size);
            merged_qual = merged_qual + qual1.substr(overlap.second + overlap_size,
                    left.size() - overlap.second - overlap_size);
        }
        else {
            merged_seq = merged_seq + rc_right.substr(overlap_size,
                    rc_right.size() - overlap_size);
            merged_qual = merged_qual + qual2.substr(overlap_size,
                    rc_right.size() - overlap_size );
        }
        VERIFY(merged_seq.size() == merged_qual.size());
        return make_pair(merged_seq, merged_qual);
    }

    std::string MergeNames(size_t index, std::string name1, std::string name2) {
        name2 = ""; // Prevent warning
        if(name1[0] == '@')
            name1 = name1.substr(1, name1.size() - 1);
        std::stringstream ss;
        ss << "@" << index << "_merged_read_" << delete_whitespaces(name1);
        return ss.str();
    }

public:
    SequenceMerger(merger_setting setting) :
        setting_(setting) { }

    FastqRead Merge(size_t index, PairedFastqRead paired_read) {
        std::string merged_name = MergeNames(index, paired_read.left_read.name,
                paired_read.right_read.name);
        std::pair<std::string, std::string> merged_seq_qual = MergeQualifiedSeq(paired_read);
        return FastqRead(merged_name, merged_seq_qual.first,
                merged_seq_qual.second);
    }
};

class PairedReadsMerger {
    SequenceMerger seq_merger_;
public:
    PairedReadsMerger(merger_setting settings) :
        seq_merger_(settings) { }

    std::vector<FastqRead> Merge(std::vector<PairedFastqRead> reads) {
        std::vector<FastqRead> merged_reads;
        double processed_rate = 10;
        for(size_t i = 0; i < reads.size(); i++) {
            FastqRead merged_read = seq_merger_.Merge(i, reads[i]);
            if(!merged_read.is_empty())
                merged_reads.push_back(merged_read);
            if(static_cast<double>(i) / static_cast<double>(reads.size()) * 100. > processed_rate) {
                INFO(processed_rate << "% reads were processed");
                processed_rate += 10;
            }
        }
        INFO("100% reads were processed");
        return merged_reads;
    }
};
