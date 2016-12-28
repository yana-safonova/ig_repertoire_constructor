#pragma once

#include "lymphocyte_type.hpp"

namespace germline_utils {

//    enum IgChainType {
//        UnknownIgChain, HeavyIgChain, KappaIgChain, LambdaIgChain
//    };
//
//    std::ostream &operator<<(std::ostream &out, const IgChainType &chain_type);
//
//    enum TcrChainType {
//        UnknownTcrChain, AlphaTcrChain, BetaTcrChain, GammaTcrChain, DeltaTcrChain
//    };
//
//    std::ostream &operator<<(std::ostream &out, const TcrChainType &tcr_chain_type);

    enum ImmuneChainType {
        UnknownImmuneChain, HeavyIgChain, KappaIgChain, LambdaIgChain,
        AlphaTcrChain, BetaTcrChain, GammaTcrChain, DeltaTcrChain
    };

    class ChainType {
        LymphocyteType lymphocyte_type_;
        ImmuneChainType chain_type_;

        void Initialize(ImmuneChainType chain_type);

    public:
        ChainType() :
                lymphocyte_type_(LymphocyteType::UnknownLymphocyte),
                chain_type_(ImmuneChainType::UnknownImmuneChain) { }

        ChainType(ImmuneChainType chain_type);

        ChainType(std::string chain_str);

        LymphocyteType Lymphocyte() const { return lymphocyte_type_; }

        ImmuneChainType Chain() const { return chain_type_; }

        bool operator==(const ChainType &obj) const {
            return lymphocyte_type_ == obj.Lymphocyte() and
                   chain_type_ == obj.Chain();
        }

        bool IsVDJ() const { return chain_type_ == ImmuneChainType::HeavyIgChain or
                    chain_type_ == ImmuneChainType::BetaTcrChain or
                    chain_type_ == ImmuneChainType::DeltaTcrChain; }
    };

    std::ostream &operator<<(std::ostream &out, const ChainType &chain_type);

    struct ChainTypeHasher {
        size_t operator()(const ChainType& chain_type) const;
    };

}