#include "logger/logger.hpp"
#include "v_recombination_event_generator.hpp"

using namespace vdj_labeler;

recombination_utils::CleavedIgGeneAlignment VRecombinationEventGenerator::GenerateCleavageEvent(
        alignment_utils::ImmuneGeneReadAlignmentPtr v_alignment,
        size_t cleavage_length)
{
    return recombination_utils::CleavedIgGeneAlignment(v_alignment,
                                                       0, // left cleavage left
                                                       static_cast<int>(cleavage_length), // right cleavage length
                                                       0, // number of SHMs in left cleavage
                                                       // number of SHMs in right cleavage
                                                       shms_calculator_.ComputeNumberSHMsForRightEvent(
                                                           v_alignment,
                                                           static_cast<int>(cleavage_length)));
}

void VRecombinationEventGenerator::GenerateCleavageEvents(
        alignment_utils::ImmuneGeneReadAlignmentPtr v_alignment,
        recombination_utils::IgGeneRecombinationEventStoragePtr v_events)
{
    // if alignment is empty cleavage can not be observed
    if(v_alignment->Empty())
        return;
    size_t min_cleavage = (v_alignment->subject().length() - 1) - v_alignment->EndSubjectPosition();
    size_t max_cleavage = std::min<size_t>(max_cleavage_, v_alignment->SubjectAlignmentLength());
    TRACE("Min cleavage: " <<  min_cleavage << ", max cleavage length in V: " << max_cleavage);
    for(size_t clen = min_cleavage; clen <= max_cleavage; clen++)
        v_events->AddEvent(GenerateCleavageEvent(v_alignment, clen));
}

recombination_utils::CleavedIgGeneAlignment VRecombinationEventGenerator::GeneratePalindromicEvent(
        alignment_utils::ImmuneGeneReadAlignmentPtr v_alignment,
        size_t palindrome_length)
{
    int event_length = int(palindrome_length) * - 1;
    return recombination_utils::CleavedIgGeneAlignment(v_alignment,
                                                       0,
                                                       event_length,
                                                       0,
                                                       // number of SHMs in right cleavage
                                                       shms_calculator_.ComputeNumberSHMsForRightEvent(v_alignment,
                                                                                                       event_length));
}

void VRecombinationEventGenerator::GeneratePalindromicEvents(
        alignment_utils::ImmuneGeneReadAlignmentPtr v_alignment,
        recombination_utils::IgGeneRecombinationEventStoragePtr v_events)
{
    if(v_alignment->EndSubjectPosition() != v_alignment->subject().length() - 1)
        return;
    size_t max_palindrome_length =
        std::min<size_t>(std::min<size_t>(max_palindrome_, v_alignment->subject().length()),
                         v_alignment->query().length() - v_alignment->EndQueryPosition() - 1);
    TRACE("Max palindrome length: " << max_palindrome_length);
    for(size_t plen = 1; plen <= max_palindrome_length; plen++)
        v_events->AddEvent(GeneratePalindromicEvent(v_alignment, plen));
}

recombination_utils::IgGeneRecombinationEventStoragePtr VRecombinationEventGenerator::ComputeEvents(
        alignment_utils::ImmuneGeneReadAlignmentPtr v_alignment)
{
    recombination_utils::IgGeneRecombinationEventStoragePtr v_events(
        new recombination_utils::IgGeneRecombinationEventStorage(germline_utils::SegmentType::VariableSegment));
    TRACE(v_alignment->Alignment());
    // generation of cleavage events
    TRACE("Generation of cleavage events");
    GenerateCleavageEvents(v_alignment, v_events);
    // generation of zero event
    TRACE("Generation of zero event");
    v_events->AddEvent(recombination_utils::CleavedIgGeneAlignment(v_alignment, 0, 0, 0, 0));
    // generation of palindromic events
    TRACE("Generartion of palindromic events");
    GeneratePalindromicEvents(v_alignment, v_events);
    return v_events;
}