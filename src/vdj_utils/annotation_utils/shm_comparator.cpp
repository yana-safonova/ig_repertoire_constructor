#include "shm_comparator.hpp"

namespace annotation_utils {
    bool SHMComparator::SHMsAreEqual(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2) {
        if(shms1.size() != shms2.size())
            return false;
        for(size_t i = 0; i < shms1.size(); i++)
            if(shms1[i] != shms2[i])
                return false;
        return true;
    }

    bool SHMComparator::SHMs1AreNestedInSHMs2(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2) {
        size_t index2 = 0;
        for(auto it1 = shms1.cbegin(); it1 != shms1.cend(); it1++) {
            bool shm_found = false;
            for(size_t i = index2; i < shms2.size(); i++)
                if(*it1 == shms2[i]) {
                    index2 = i + 1;
                    shm_found = true;
                    break;
                }
            if(!shm_found)
                return false;
        }
        return true;
    }

    size_t SHMComparator::GetNumberOfIntersections(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2) {
        size_t num_shared_shms = 0;
        size_t index2 = 0;
        for(auto it1 = shms1.cbegin(); it1 != shms1.cend(); it1++) {
            for(size_t i = index2; i < shms2.size(); i++)
                if(*it1 == shms2[i]) {
                    num_shared_shms++;
                    index2 = i + 1;
                    break;
                }
        }
        return num_shared_shms;
    }

    bool SHMComparator::AddedSHMsAreSynonimous(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2) {
        size_t index1 = 0;
        for(auto it2 = shms2.cbegin(); it2 != shms2.cend(); it2++) {
            bool shm_found = false;
            for(size_t i = index1; i < shms1.size(); i++)
                if(*it2 == shms1[i]) {
                    index1 = i + 1;
                    shm_found = true;
                    break;
                }
            if(!shm_found && !it2->IsSynonymous()) {
                return false;
            }
        }
        return true;
    }
}
