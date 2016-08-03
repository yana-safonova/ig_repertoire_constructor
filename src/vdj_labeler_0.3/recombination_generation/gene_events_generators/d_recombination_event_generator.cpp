#include "logger/logger.hpp"
#include "d_recombination_event_generator.hpp"

using namespace vdj_labeler;

int DRecombinationEventGenerator::ComputeMinLeftBound(
        const alignment_utils::ImmuneGeneReadAlignment &d_alignment) const
{
    if(d_alignment.StartSubjectPosition() != 0)
        return int(d_alignment.StartSubjectPosition());
    return int(max_palindrome_) * -1;
}

int DRecombinationEventGenerator::ComputeMinRightBound(const alignment_utils::ImmuneGeneReadAlignment &d_alignment) const
{
    // if cleavage is occurred at the right end of D segment, min right bound is cleavage size from alignment
    // TRACE("if(d_alignment.EndSubjectPosition() != d_alignment.Subject().length() - 1)" << d_alignment.EndSubjectPosition() << d_alignment.Subject().length());
    if(d_alignment.EndSubjectPosition() != d_alignment.Subject().length() - 1)
        return int(d_alignment.Subject().length() - d_alignment.EndSubjectPosition() - 1);
    // if gene segment is not cleaved, min left bound is max palindrome length
    return int(max_palindrome_) * -1;
}

size_t DRecombinationEventGenerator::ComputeMaxRightConsistentCleavage(
        const alignment_utils::ImmuneGeneReadAlignment &d_alignment,
        const int left_event_size) const
{
    size_t min_right_cleavage = d_alignment.Subject().length() - d_alignment.EndSubjectPosition() - 1;
    size_t read_cleavage = 0;
    if(left_event_size > 0) {
        assert(size_t(left_event_size) >= d_alignment.StartSubjectPosition());
        read_cleavage = left_event_size - d_alignment.StartSubjectPosition();
    }
    return size_t(std::min(static_cast<int>(max_cleavage_),
                           static_cast<int>(d_alignment.QueryAlignmentLength())) -
                                    read_cleavage + min_right_cleavage);
}

void DRecombinationEventGenerator::GenerateRightConsistentEvents(
        const alignment_utils::ImmuneGeneReadAlignment &d_alignment,
        const int left_event_size,
        recombination_utils::IgGeneRecombinationEventStorage &d_events) const
{
    int min_right_bound = ComputeMinRightBound(d_alignment);
    int max_right_bound = static_cast<int>(ComputeMaxRightConsistentCleavage(d_alignment, left_event_size));
    TRACE("Left bound of right events: " << min_right_bound);
    TRACE("Right bound of right events: " << max_right_bound);
    for(int relen = min_right_bound; relen <= max_right_bound; relen++) {
        TRACE("== Right current event: " << relen << ".");
        d_events.AddEvent(recombination_utils::CleavedIgGeneAlignment(d_alignment,
                                                                      left_event_size,
                                                                      relen,
                                                                      shm_calculator_.ComputeNumberSHMsForLeftEvent(
                                                                          d_alignment, left_event_size),
                                                                      shm_calculator_.ComputeNumberSHMsForRightEvent(
                                                                          d_alignment, relen)));
    }
}

recombination_utils::IgGeneRecombinationEventStorage DRecombinationEventGenerator::ComputeEvents(
        const alignment_utils::ImmuneGeneReadAlignment &d_alignment) const
{
    recombination_utils::IgGeneRecombinationEventStorage d_events(germline_utils::SegmentType::DiversitySegment);
    if(d_alignment.Empty())
        return d_events;
    int min_left_bound = ComputeMinLeftBound(d_alignment);
    int max_left_bound =
        static_cast<int>(std::min<size_t>(d_alignment.AlignmentLength() + d_alignment.StartSubjectPosition(),
                                          max_cleavage_));
    INFO(d_alignment.Alignment());
    INFO(d_alignment.Subject());
    INFO(d_alignment.Query());
    INFO("Left bound of left events: " << min_left_bound);
    INFO("Right bound of left events: " << max_left_bound);
    //  we iterate from max allowed palindrome to max allowed cleavage and
    // consider that this event occurred at the start of D segment
    for(int elen = min_left_bound; elen <= max_left_bound; elen++) {
        TRACE("==== Left current event: " << elen << "...");
        GenerateRightConsistentEvents(d_alignment, elen, d_events);
    }
    return d_events;
}
