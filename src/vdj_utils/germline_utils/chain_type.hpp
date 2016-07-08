#pragma once

#include "lymphocyte_type.hpp"

namespace germline_utils {

    enum IgChainType {
        UnknownIgChain, HeavyIgChain, KappaIgChain, LambdaIgChain
    };

    std::ostream &operator<<(std::ostream &out, const IgChainType &chain_type);

    enum TcrChainType {
        UnknownTcrChain, AlphaTcrChain, BetaTcrChain, GammaTcrChain, DeltaTcrChain
    };

    std::ostream &operator<<(std::ostream &out, const TcrChainType &tcr_chain_type);

    enum ImmuneChainType {
        UnknownImmuneChain, HeavyIgChain, KappaIgChain, LambdaIgChain,
        AlphaTcrChain, BetaTcrChain, GammaTcrChain, DeltaTcrChain
    };

    class ChainType {
        LymphocyteType lymphocyte_type_;
        ImmuneChainType chain_type_;

    public:
        ChainType() :
                lymphocyte_type_(LymphocyteType::UnknownLymphocyte),
                chain_type_(chain_type_::UnknownChain) { }

        ChainType(ImmuneChainType chain_type);

        LymphocyteType Lymphocyte() const { return lymphocyte_type_; }

        ImmuneChainType Chain() const { return chain_type_; }

        bool operator==(const ChainType &obj) {
            return lymphocyte_type_ == obj.Lymphocyte() and
                   chain_type_ == obj.Chain();
        }
    };

    std::ostream &operator<<(std::ostream &out, ChainType &chain_type);

}