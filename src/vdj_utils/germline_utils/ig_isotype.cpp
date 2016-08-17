#include "ig_isotype.hpp"

IgIsotype::InnerIgIsotype IgIsotype::ConvertStringToIgIsotype(std::string isotype_str) {
    if(isotype_str == "IGG")
        return InnerIgIsotype::IgG;
    if(isotype_str == "IGM")
        return InnerIgIsotype::IgM;
    if(isotype_str == "IGA")
        return InnerIgIsotype::IgA;
    if(isotype_str == "IGD")
        return InnerIgIsotype::IgD;
    if(isotype_str == "IGE")
        return InnerIgIsotype::IgE;
    if(isotype_str == "IGK")
        return InnerIgIsotype::IgK;
    if(isotype_str == "IGL")
        return InnerIgIsotype::IgL;
    return InnerIgIsotype::unknown_isotype;
}