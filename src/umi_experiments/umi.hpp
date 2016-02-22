#include <seqan/seq_io.h>
#include <verify.hpp>
#include "../ig_tools/utils/string_tools.hpp"

class Umi {
public:
    explicit Umi(const seqan::Dna5String& umi) : umi_(umi) {}

    bool operator==(const Umi &other) const { return umi_ == other.umi_; }

    seqan::Dna5String GetString() const { return umi_; }

private:
    const seqan::Dna5String umi_;
};

namespace std {
    template<>
    struct hash<seqan::Dna5String> {
        size_t operator()(const seqan::Dna5String& str) const {
            size_t h = 0;
            for (auto c : str) {
                h = h * 31 + seqan::ordValue(c);
            }
            return h;
        }
    };

    template<>
    struct hash<Umi> {
        size_t operator()(const Umi& umi) const {
            size_t h = hash<seqan::Dna5String>()(umi.GetString());
            return h;
        }
    };
}

// TODO: replace its usage with reading another file with umi's: duplicates umi_to_fastq
void extract_umi(const std::vector<seqan::CharString>& ids, std::vector<seqan::Dna5String>& umis, std::vector<seqan::DnaQString>& umi_quals) {
    for (auto& id : ids) {
        auto id_string = seqan::toCString(id);
        auto split_by_umi = split(id_string, "UMI");
        VERIFY_MSG(split_by_umi.size() == 2, "Either no UMI info, or too much of it in meta: " << id_string);
        auto umi_info = split_by_umi[1].substr(1);
        auto colon = umi_info.find(":");
        VERIFY_MSG(colon != std::string::npos, "Can't parse UMI info: " << umi_info);
        auto umi = umi_info.substr(0, colon);
        auto qual = umi_info.substr(colon + 1);
        VERIFY_MSG(umi.length() == qual.length(), "UMI and its quality are of different lengths: " << umi_info);
        umis.push_back(umi);
        umi_quals.push_back(qual);
    }
}
