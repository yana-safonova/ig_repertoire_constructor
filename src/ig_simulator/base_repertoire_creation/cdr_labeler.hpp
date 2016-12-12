#pragma once

#include "../ig_structs/cdr_settings.hpp"
#include "../repertoire.hpp"

template<class IgVariableRegionPtr, class CDRLabelingStrategy, class CDRSettings>
class CDRLabeler {
    CDRLabelingStrategy labeling_strategy_;
public:
    CDRLabeler(CDRLabelingStrategy labeling_strategy) :
            labeling_strategy_(labeling_strategy) { }

    IgVariableRegionPtr LabelCDRs(IgVariableRegionPtr ig_variable_region_ptr) {
        CDRSettings cdr_settings = labeling_strategy_.LabelCDRs(ig_variable_region_ptr);
        ig_variable_region_ptr->SetCDRSettings(cdr_settings);
        return ig_variable_region_ptr;
    }
};

// ----------------------------------------------------------
// 30 to 36, 49 to 65, and 95 to 103. Total: 120
struct HC_CDRConstants {
    static constexpr double cdr1_start = .25;
    static constexpr double cdr1_end = .3;
    static constexpr double cdr2_start = .41;
    static constexpr double cdr2_end = .54;
    static constexpr double cdr3_start = .79;
    static constexpr double cdr3_end = .86;
};

// 28 to 35, 49 to 59, and 92 to 103. Total: 120
struct LC_CDRConstants {
    static constexpr double cdr1_start = .23;
    static constexpr double cdr1_end = .29;
    static constexpr double cdr2_start = .41;
    static constexpr double cdr2_end = .49;
    static constexpr double cdr3_start = .77;
    static constexpr double cdr3_end = .86;
};

template<class IgVariableRegionPtr, class CDRSettings>
class HC_SimpleCDRLabelingStrategy {
public:
    CDRSettings LabelCDRs(IgVariableRegionPtr ig_variable_region_ptr) {
        size_t vregion_length = ig_variable_region_ptr->Length();
        CDRLabel cdr1(size_t(double(vregion_length) * HC_CDRConstants::cdr1_start),
                      size_t(double(vregion_length) * HC_CDRConstants::cdr1_end));
        CDRLabel cdr2(size_t(double(vregion_length) * HC_CDRConstants::cdr2_start),
                      size_t(double(vregion_length) * HC_CDRConstants::cdr2_end));
        CDRLabel cdr3(size_t(double(vregion_length) * HC_CDRConstants::cdr3_start),
                      size_t(double(vregion_length) * HC_CDRConstants::cdr3_end));
        return CDRSettings(cdr1, cdr2, cdr3);
    }
};

template<class IgVariableRegionPtr, class CDRSettings>
class LC_SimpleCDRLabelingStrategy {
public:
    CDRSettings LabelCDRs(IgVariableRegionPtr ig_variable_region_ptr) {
        size_t vregion_length = ig_variable_region_ptr->Length();
        CDRLabel cdr1(size_t(double(vregion_length) * LC_CDRConstants::cdr1_start),
                      size_t(double(vregion_length) * LC_CDRConstants::cdr1_end));
        CDRLabel cdr2(size_t(double(vregion_length) * LC_CDRConstants::cdr2_start),
                      size_t(double(vregion_length) * LC_CDRConstants::cdr2_end));
        CDRLabel cdr3(size_t(double(vregion_length) * LC_CDRConstants::cdr3_start),
                      size_t(double(vregion_length) * LC_CDRConstants::cdr3_end));
        return CDRSettings(cdr1, cdr2, cdr3);
    }
};

typedef HC_SimpleCDRLabelingStrategy<HC_VariableRegionPtr, CDRSettings> HC_CDRLabelingStrategy;
typedef LC_SimpleCDRLabelingStrategy<LC_VariableRegionPtr, CDRSettings> LC_CDRLabelingStrategy;

// ----------------------------------------------------------
typedef CDRLabeler<HC_VariableRegionPtr, HC_CDRLabelingStrategy, CDRSettings> HC_CDRLabeler;
typedef CDRLabeler<LC_VariableRegionPtr, LC_CDRLabelingStrategy, CDRSettings> LC_CDRLabeler;