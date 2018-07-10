#pragma once

#include "shm_annotation/shm_annotation.hpp"

namespace annotation_utils {
    class SHMComparator {
    public:
        static bool SHMsAreEqual(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2);

        static bool SHMs1AreNestedInSHMs2(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2);

        static bool AllSHMs1InsertionBlocksArePresentedInSHMs2(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2);

        static bool SHMsInsertionBlocksAreEqual(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2);

        static size_t GetNumberOfIntersections(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2);

        static bool AddedSHMsAreSynonimous(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2);

        static bool AllAddedSHMs1HaveIdenticallyPositionedSHMs2(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2);

        static bool IndividualSHMsAreIdenticallyPositioned(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2);

        static std::vector<SHM> GetAddedSHMs(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2);
    };
}
