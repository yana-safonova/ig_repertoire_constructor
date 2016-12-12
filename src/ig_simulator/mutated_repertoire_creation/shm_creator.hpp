#pragma once

#include "../ig_structs/ig_structure_structs.hpp"
#include "../repertoire.hpp"
#include "shm_strategies/rgyw_wrcy_strategy.hpp"
#include "shm_strategies/cdr_based_random_strategy.hpp"

template<class IgVariableRegionPtr, class SHMCreationStrategy, class SHMCreationSettings>
class SHMCreator {
    SHMCreationStrategy strategy_;
public:
    SHMCreator(SHMCreationStrategy strategy) :
        strategy_(strategy) { }

    IgVariableRegionPtr CreateSHM(IgVariableRegionPtr ig_variable_region_ptr) {
        SHMCreationSettings settings =  strategy_.CreateSHM(ig_variable_region_ptr);
        ig_variable_region_ptr->SetSHMSettings(settings);
        return ig_variable_region_ptr;
    }
};

// ----------------------------------------------------------
// composite strategy
template<class IgVariableRegionPtr, class SHMSettings>
class CompositeSHMCreationStrategy {
    RgywWrcySHMStrategy<IgVariableRegionPtr, SHMSettings> rgyw_wrcy_strategy_;
    CDRBasedRandomSHMStrategy<IgVariableRegionPtr, SHMSettings> cdr_based_strategy_;
public:
    CompositeSHMCreationStrategy(RgywWrcySHMStrategy<IgVariableRegionPtr, SHMSettings> rgyw_wrcy_strategy,
                                 CDRBasedRandomSHMStrategy<IgVariableRegionPtr, SHMSettings> cdr_based_strategy) :
        rgyw_wrcy_strategy_(rgyw_wrcy_strategy),
        cdr_based_strategy_(cdr_based_strategy) { }

    SHMSettings CreateSHM(IgVariableRegionPtr ig_variable_region_ptr) {
        SHMSettings shm_settings1 = rgyw_wrcy_strategy_.CreateSHM(ig_variable_region_ptr);
        SHMSettings shm_settings2 = cdr_based_strategy_.CreateSHM(ig_variable_region_ptr);
        for(auto it = shm_settings2.begin(); it != shm_settings2.end(); it++)
            shm_settings1.Add(it->second);
        return shm_settings1;
    }
};
// ---------------------------------------------------------

typedef CompositeSHMCreationStrategy<HC_VariableRegionPtr, SHMSettings> HC_CompositeSHMCreationStrategy;
typedef CompositeSHMCreationStrategy<LC_VariableRegionPtr, SHMSettings> LC_CompositeSHMCreationStrategy;

typedef SHMCreator<HC_VariableRegionPtr, HC_CompositeSHMCreationStrategy, SHMSettings> HC_SHMCreator;
typedef SHMCreator<LC_VariableRegionPtr, LC_CompositeSHMCreationStrategy, SHMSettings> LC_SHMCreator;