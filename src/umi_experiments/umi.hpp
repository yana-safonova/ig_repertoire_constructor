#include <seqan/seq_io.h>

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
