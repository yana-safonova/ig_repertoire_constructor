#pragma once

#include "vdj_alignments/vdj_hits.hpp"

#include "recombination_generator.hpp"
#include "recombination_utils/recombination.hpp"
#include "recombination_utils/recombination_storage.hpp"
#include "gene_events_generators/ig_gene_recombination_event_generator.hpp"
#include "insertion_events_generators/insertion_event_generator.hpp"

namespace vdj_labeler {

class CustomHeavyChainRecombinationGenerator :
        public AbstractRecombinationGenerator<recombination_utils::HCRecombination,
                                              recombination_utils::HcRecombinationStorage>
{
    IgGeneRecombinationEventsGenerator &v_events_generator_;
    IgGeneRecombinationEventsGenerator &d_events_generator_;
    IgGeneRecombinationEventsGenerator &j_events_generator_;
    InsertionEventGenerator& vd_insertion_generator_;
    InsertionEventGenerator& dj_insertion_generator_;

    // inner structures
    std::vector<recombination_utils::IgGeneRecombinationEventStorage> v_storages_;
    std::vector<recombination_utils::IgGeneRecombinationEventStorage> d_storages_;
    std::vector<recombination_utils::IgGeneRecombinationEventStorage> j_storages_;

    void Clear();

    void ComputeVEventStorages(const VDJHits &vdj_hits);

    void ComputeDEventStorages(const VDJHits &vdj_hits);

    void ComputeJEventStorages(const VDJHits &vdj_hits);

    // if d alignment is empty, we will assign length from V end to J start to the first insertion
    // and zero length to the second insertion
    // todo: extract to separate class
    std::pair<recombination_utils::NongenomicInsertion, recombination_utils::NongenomicInsertion>
            RefineNongenomicInsertions(const recombination_utils::NongenomicInsertion &vd_insertion,
                                       const recombination_utils::NongenomicInsertion &dj_insertion) const;

    void CreateRecombinations(
        recombination_utils::HcRecombinationStorage &recombination_storage,
        const recombination_utils::CleavedIgGeneAlignment &v_gene,
        const recombination_utils::CleavedIgGeneAlignment &d_gene,
        const recombination_utils::CleavedIgGeneAlignment &j_gene,
        recombination_utils::InsertionEventStorage &vd_insertions,
        recombination_utils::InsertionEventStorage &dj_insertions) const;

    void CreateRecombinations(
        recombination_utils::HcRecombinationStorage &recombination_storage,
        const recombination_utils::IgGeneRecombinationEventStorage &v_events,
        const recombination_utils::IgGeneRecombinationEventStorage &d_events,
        const recombination_utils::IgGeneRecombinationEventStorage &j_events) const;

public:
    CustomHeavyChainRecombinationGenerator(IgGeneRecombinationEventsGenerator &v_events_generator,
                                           IgGeneRecombinationEventsGenerator &d_events_generator,
                                           IgGeneRecombinationEventsGenerator &j_events_generator,
                                           InsertionEventGenerator& vd_insertion_generator,
                                           InsertionEventGenerator& dj_insertion_generator) :
            v_events_generator_(v_events_generator),
            d_events_generator_(d_events_generator),
            j_events_generator_(j_events_generator),
            vd_insertion_generator_(vd_insertion_generator),
            dj_insertion_generator_(dj_insertion_generator) { }

    recombination_utils::HcRecombinationStorage ComputeRecombinations(const VDJHits &vdj_hits) override;
};

} // End namespace vdj_labeler