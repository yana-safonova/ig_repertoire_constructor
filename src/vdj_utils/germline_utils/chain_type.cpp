#include <verify.hpp>

#include "chain_type.hpp"

#include <vector>
#include <ostream>

namespace germline_utils {

//    std::ostream &operator<<(std::ostream &out, const IgChainType &ig_chain_type) {
//        if (ig_chain_type == IgChainType::HeavyIgChain)
//            out << "IGH";
//        else if (ig_chain_type == IgChainType::KappaIgChain)
//            out << "IGK";
//        else if (ig_chain_type == IgChainType::LambdaIgChain)
//            out << "IGL";
//        else
//            out << "Unknown Ig chain";
//        return out;
//    }
//
//    std::ostream &operator<<(std::ostream &out, const TcrChainType &tcr_chain_type) {
//        if (tcr_chain_type == TcrChainType::AlphaTcrChain)
//            out << "TRA";
//        else if (tcr_chain_type == TcrChainType::BetaTcrChain)
//            out << "TRB";
//        else if (tcr_chain_type == TcrChainType::GammaTcrChain)
//            out << "TRG";
//        else if (tcr_chain_type == TcrChainType::DeltaTcrChain)
//            out << "TRD";
//        else
//            out << "Unknown TCR chain";
//        return out;
//    }

    std::vector<std::string> immune_chain_strs { "Unknown chain type", "IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD" };

    void check_chain_str_correctness_fatal(std::string chain_str) {
        bool chain_str_was_not_found = false;
        for(auto it = immune_chain_strs.begin(); it != immune_chain_strs.end(); it++)
            if(*it == chain_str)
                return;
        VERIFY_MSG(chain_str_was_not_found, "Chain type " << chain_str << " was not recognized");
    }

    ImmuneChainType get_chain_type_by_str(std::string chain_str) {
        for(size_t i = 0; i < immune_chain_strs.size(); i++)
            if(immune_chain_strs[i] == chain_str)
                return ImmuneChainType(i);
        return ImmuneChainType::UnknownImmuneChain;
    }

    std::ostream &operator<<(std::ostream &out, const ImmuneChainType &chain_type) {
        out << immune_chain_strs[int(chain_type)];
        return out;
    }

    bool immune_chain_is_tcr(ImmuneChainType immune_chain_type) {
        return immune_chain_type == ImmuneChainType::AlphaTcrChain or
               immune_chain_type == ImmuneChainType::BetaTcrChain or
               immune_chain_type == ImmuneChainType::GammaTcrChain or
               immune_chain_type == ImmuneChainType::DeltaTcrChain;
    }

    bool immune_chain_is_ig(ImmuneChainType immune_chain_type) {
        return immune_chain_type == ImmuneChainType::HeavyIgChain or
               immune_chain_type == ImmuneChainType::KappaIgChain or
               immune_chain_type == ImmuneChainType::LambdaIgChain;
    }

    void ChainType::Initialize(ImmuneChainType chain_type) {
        if (immune_chain_is_tcr(chain_type))
            lymphocyte_type_ = LymphocyteType::TLymphocyte;
        else if (immune_chain_is_ig(chain_type))
            lymphocyte_type_ = LymphocyteType::BLymphocyte;
        else
            lymphocyte_type_ = LymphocyteType::UnknownLymphocyte;
        chain_type_ = chain_type;
    }

    ChainType::ChainType(ImmuneChainType chain_type) {
        Initialize(chain_type);
    }

    ChainType::ChainType(std::string chain_str) {
        check_chain_str_correctness_fatal(chain_str);
        ImmuneChainType chain_type = get_chain_type_by_str(chain_str);
        Initialize(chain_type);
    }

    std::ostream &operator<<(std::ostream &out, const ChainType &chain_type) {
        out << chain_type.Chain();
        return out;
    }

    size_t ChainTypeHasher::operator()(const ChainType &chain_type) const {
        return std::hash<int>()(chain_type.Chain()) * std::hash<int>()(chain_type.Lymphocyte());
    }
}