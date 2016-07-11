#include "logger/logger.hpp"
#include "d_recombination_event_generator.hpp"

using namespace std;

int DRecombinationEventGenerator::ComputeMinLeftBound(IgGeneAlignmentPtr d_alignment) {
    if(d_alignment->Positions().GeneStartPos() != 0)
        return int(d_alignment->Positions().GeneStartPos());
    return int(max_palindrome_) * -1;
}

int DRecombinationEventGenerator::ComputeMinRightBound(IgGeneAlignmentPtr d_alignment) {
    // if cleavage is occurred at the right end of D segment, min right bound is cleavage size from alignment
    if(d_alignment->Positions().GeneEndPos() != d_alignment->GeneLength() - 1)
        return int(d_alignment->GeneLength() - d_alignment->Positions().GeneEndPos() - 1);
    // if gene segment is not cleaved, min left bound is max palindrome length
    return int(max_palindrome_) * -1;
}

size_t DRecombinationEventGenerator::ComputeMaxRightConsistentCleavage(IgGeneAlignmentPtr d_alignment,
                                                                    int left_event_size) {
    size_t min_right_cleavage = d_alignment->GeneLength() - d_alignment->Positions().GeneEndPos() - 1;
    size_t read_cleavage = 0;
    if(left_event_size > 0) {
        assert(size_t(left_event_size) >= d_alignment->Positions().GeneStartPos());
        read_cleavage = left_event_size - d_alignment->Positions().GeneStartPos();
    }
    return size_t(min<int>(int(max_cleavage_),
                           int(d_alignment->ReadAlignmentLength())) - read_cleavage + min_right_cleavage);
}

void DRecombinationEventGenerator::GenerateRightConsistentEvents(IgGeneAlignmentPtr d_alignment, int left_event_size,
                                                                 IgGeneRecombinationEventStoragePtr d_events) {
    int min_right_bound = ComputeMinRightBound(d_alignment);
    int max_right_bound = int(ComputeMaxRightConsistentCleavage(d_alignment, left_event_size));
    TRACE("Left bounf of right events: " << min_right_bound);
    TRACE("Right bound of right events: " << max_right_bound);
    for(int relen = min_right_bound; relen <= max_right_bound; relen++) {
        TRACE("== Right current event: " << relen << ".");
        d_events->AddEvent(CleavedIgGeneAlignment(d_alignment,
                                                  left_event_size,
                                                  relen,
                                                  shm_calculator_.ComputeNumberSHMsForLeftEvent(d_alignment,
                                                                                                left_event_size),
                                                  shm_calculator_.ComputeNumberSHMsForRightEvent(d_alignment, relen)));
    }
}

IgGeneRecombinationEventStoragePtr DRecombinationEventGenerator::ComputeEvents(IgGeneAlignmentPtr d_alignment) {
    IgGeneRecombinationEventStoragePtr d_events(new IgGeneRecombinationEventStorage(IgGeneType::diversity_gene));
    if(d_alignment->IsEmpty())
        return d_events;
    int min_left_bound = ComputeMinLeftBound(d_alignment);
    int max_left_bound = int(min<size_t>(d_alignment->ReadAlignmentLength() + d_alignment->Positions().GeneStartPos(),
                                         max_cleavage_));
    TRACE(*d_alignment);
    TRACE("Left bound of left events: " << min_left_bound);
    TRACE("Right bound of left events: " << max_left_bound);
    //  we iterate from max allowed palindrome to max allowed cleavage and
    // consider that this event occurred at the start of D segment
    for(int elen = min_left_bound; elen <= max_left_bound; elen++) {
        TRACE("==== Left current event: " << elen << "...");
        GenerateRightConsistentEvents(d_alignment, elen, d_events);
    }
    return d_events;
}
