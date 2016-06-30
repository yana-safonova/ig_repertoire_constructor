#include <boost/algorithm/string/trim.hpp>
#include <boost/format.hpp>
#include <unordered_map>
#include <logger/log_writers.hpp>
#include "umi_utils.hpp"

void extract_barcodes_from_read_ids(const std::vector<seqan::CharString>& input_ids, std::vector<seqan::Dna5String>& umis) {
    std::vector<seqan::DnaQString> umi_quals;
    extract_barcodes_from_read_ids(input_ids, umis, umi_quals);
}

void extract_barcodes_from_read_ids(const std::vector<seqan::CharString>& input_ids, std::vector<seqan::Dna5String>& umis,
                                    std::vector<seqan::DnaQString>& umi_quals) {
    for (auto& id : input_ids) {
        std::string s = seqan_string_to_string(id);
        auto split_by_umi = split(s, "UMI");
        VERIFY_MSG(split_by_umi.size() <= 2, "Too much 'UMI' strings in read id");
        if (split_by_umi.size() == 1) {
            split_by_umi = split(s, "BARCODE");
            VERIFY_MSG(split_by_umi.size() > 1, boost::format("Could not find both 'UMI' and 'BARCODE' in read id '%s'") % s);
            VERIFY_MSG(split_by_umi.size() < 3, "Too much 'BARCODE' strings in read id");
        }
        std::string meta = split_by_umi[0];
        boost::algorithm::trim(meta);
        std::string umi_info = split_by_umi[1].substr(1);
        VERIFY(!umi_info.empty());
        size_t colon = umi_info.find(':');
        if (colon == std::string::npos) {
            umis.push_back(umi_info);
        } else {
            auto umi = umi_info.substr(0, colon);
            auto qual = umi_info.substr(colon + 1);
            VERIFY_MSG(umi.length() == qual.length(), "UMI and its quality are of different lengths: " << umi_info);
            umis.push_back(umi);
            umi_quals.push_back(qual);
        }
        VERIFY_MSG(umi_quals.size() == umis.size() || umi_quals.size() == 0, "Found both UMIs with quality data and without.");
    }
}

void group_reads_by_umi(const std::vector<seqan::Dna5String>& umis, std::unordered_map<Umi, std::vector<size_t> >& umi_to_reads) {
    INFO("Grouping reads by barcodes");
    for (size_t i = 0; i < umis.size(); i ++) {
        Umi umi(umis[i]);
        auto& reads = umi_to_reads[umi];
        reads.push_back(i);
    }

    size_t reads_count = 0;
    for (const auto& entry : umi_to_reads) {
        reads_count += entry.second.size();
    }
    VERIFY_MSG(reads_count == umis.size(), "Expected to have " << umis.size() << " reads corresponding to some UMI, but got " << reads_count);

    INFO(umi_to_reads.size() << " unique barcodes found");
}
