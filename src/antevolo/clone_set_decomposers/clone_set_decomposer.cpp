#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include "clone_set_decomposer.hpp"

namespace  antevolo {
    std::string CloneSetDecomposer::GetGeneBaseName(seqan::CharString name) const {
        std::string gene_name = std::string(seqan::toCString(name));
        //return gene_name;
        std::vector<std::string> splits;
        boost::split(splits, gene_name, boost::is_any_of("*"), boost::token_compress_on);
        return splits[0];
    }
}