#pragma once

#include <string>



class IgIsotype {
    enum InnerIgIsotype  { unknown_isotype,
        IgG,
        IgM,
        IgA,
        IgE,
        IgD,
        IgK,
        IgL
    };

    InnerIgIsotype isotype_;
    std::string isotype_str_;

    InnerIgIsotype ConvertStringToIgIsotype(std::string isotype_str);

public:
    IgIsotype() :
            isotype_(InnerIgIsotype::unknown_isotype),
            isotype_str_() { }

    IgIsotype(std::string isotype_str) :
            isotype_(ConvertStringToIgIsotype(isotype_str)),
            isotype_str_(isotype_str) { }

    bool IsHeavyChain() const { return isotype_ == InnerIgIsotype::IgG or isotype_ == InnerIgIsotype::IgM or
                                       isotype_ == InnerIgIsotype::IgD or isotype_ == InnerIgIsotype::IgE or
                                       isotype_ == InnerIgIsotype::IgA; }

    bool IsLightChain() const { return isotype_ == InnerIgIsotype::IgK or isotype_ == InnerIgIsotype::IgL; };

    bool IsIgG() const { return isotype_ == InnerIgIsotype::IgG; }

    bool IsIgM() const { return isotype_ == InnerIgIsotype::IgM; }

    bool IsIgD() const { return isotype_ == InnerIgIsotype::IgD; }

    bool IsIgE() const { return isotype_ == InnerIgIsotype::IgE; }

    bool IsIgA() const { return isotype_ == InnerIgIsotype::IgA; }

    bool IsIgL() const { return isotype_ == InnerIgIsotype::IgL; }

    bool IsIgK() const { return isotype_ == InnerIgIsotype::IgK; }

    bool IsUnknown() const { return isotype_ == InnerIgIsotype::unknown_isotype; }

    bool operator==(const IgIsotype &ig_isotype) const { return ig_isotype.isotype_ == isotype_; }

    std::string str() const { return isotype_str_; }

    bool operator<(const IgIsotype &ig_isotype) const { return isotype_ < ig_isotype.isotype_; }
};

class IgIsotypeHelper {
public:
    static IgIsotype GetKappaIsotype() { return IgIsotype("IGK"); }

    static IgIsotype GetLambdaIsotype() { return IgIsotype("IGL"); }
};