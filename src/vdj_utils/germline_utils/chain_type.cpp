#include "chain_type.hpp"

namespace germline_utils {

    std::ostream &operator<<(std::ostream &out, const IgChainType &ig_chain_type) {
        if (ig_chain_type == IgChainType::HeavyIgChain)
            out << "IGH";
        else if (ig_chain_type == IgChainType::KappaIgChain)
            out << "IGK";
        else if (ig_chain_type == IgChainType::LambdaIgChain)
            out << "IGL";
        else
            out << "Unknown Ig chain";
        return out;
    }

    std::ostream &operator<<(std::ostream &out, const TcrChainType &tcr_chain_type) {
        if (tcr_chain_type == TcrChainType::AlphaTcrChain)
            out << "TRA";
        else if (tcr_chain_type == TcrChainType::BetaTcrChain)
            out << "TRB";
        else if (tcr_chain_type == TcrChainType::GammaTcrChain)
            out << "TRG";
        else if (tcr_chain_type == TcrChainType::DeltaTcrChain)
            out << "TRD";
        else
            out << "Unknown TCR chain";
        return out;
    }

    std::ostream &operator<<(std::ostream &out, ImmuneChainType &chain_type) {
        if (chain_type == ImmuneChainType::AlphaTcrChain)
            out << "TRA";
        else if (chain_type == ImmuneChainType::BetaTcrChain)
            out << "TRB";
        else if (chain_type == ImmuneChainType::GammaTcrChain)
            out << "TRG";
        else if (chain_type == ImmuneChainType::DeltaTcrChain)
            out << "TRD";
        else if (chain_type == ImmuneChainType::HeavyIgChain)
            out << "IGH";
        else if (chain_type == ImmuneChainType::KappaIgChain)
            out << "IGK";
        else if (chain_type == ImmuneChainType::LambdaIgChain)
            out << "IGL";
        else
            out << "Unknown chain type";
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

    ChainType::ChainType(ImmuneChainType chain_type) {
        if (immune_chain_is_tcr(chain_type))
            lymphocype_type_ = LymphocyteType::TLymphocyte;
        else if (immune_chain_is_ig(chain_type))
            lymphocype_type_ = LymphocyteType::BLymphocyte;
        lymphocype_type_ = LymphocyteType
    }

    std::ostream &operator<<(std::ostream &out, ChainType &chain_type) {
        out << chain_type.ChainType();
        return out;
    }

}