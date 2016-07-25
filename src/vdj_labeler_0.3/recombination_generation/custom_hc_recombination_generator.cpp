#include "logger/logger.hpp"
#include "custom_hc_recombination_generator.hpp"

#include <iostream>

using namespace vdj_labeler;

void CustomHeavyChainRecombinationGenerator::Clear() {
    v_storages_.clear();
    d_storages_.clear();
    j_storages_.clear();
}

void CustomHeavyChainRecombinationGenerator::ComputeVEventStorages(const VDJHits &vdj_hits) {
    TRACE("Computation of event vector for V hits starts...");
    size_t v_events_num = 0;
    for(size_t vi = 0; vi < vdj_hits.VHitsNumber(); vi++) {
        auto &v_alignment = vdj_hits.GetAlignmentByIndex(germline_utils::SegmentType::VariableSegment, vi);
        auto v_events = v_events_generator_.ComputeEvents(v_alignment);
        v_events_num += v_events.size();
        v_storages_.emplace_back(std::move(v_events));
    }
    TRACE(v_events_num << " events were computed for " << vdj_hits.VHitsNumber() << " V hits");
}

void CustomHeavyChainRecombinationGenerator::ComputeDEventStorages(const VDJHits &vdj_hits) {
    TRACE("Computation of events vector for D hits starts...");
    size_t d_events_num = 0;
    for(size_t di = 0; di < vdj_hits.DHitsNumber(); di++) {
        auto &d_alignment = vdj_hits.GetAlignmentByIndex(germline_utils::SegmentType::DiversitySegment, di);
        auto d_events = d_events_generator_.ComputeEvents(d_alignment);
        d_events_num += d_events.size();
        d_storages_.emplace_back(std::move(d_events));
    }
    TRACE(d_events_num << " events were computed for " << vdj_hits.DHitsNumber() << " D hits");
    if(d_events_num == 0) {
        assert(vdj_hits.DHitsNumber() == 1);
        d_storages_[0].AddEvent(recombination_utils::CleavedIgGeneAlignment(
                vdj_hits.GetAlignmentByIndex(germline_utils::SegmentType::DiversitySegment, 0), 0, 0, 0, 0));
        TRACE("Himeric event was added to D event storage");
    }
}

void CustomHeavyChainRecombinationGenerator::ComputeJEventStorages(const VDJHits &vdj_hits) {
    size_t j_events_num = 0;
    TRACE("Computation of events vector for J segments starts...");
    for(size_t ji = 0; ji < vdj_hits.JHitsNumber(); ji++) {
        TRACE(ji);
        auto &j_alignment = vdj_hits.GetAlignmentByIndex(germline_utils::SegmentType::JoinSegment, ji);
        TRACE(ji);
        auto j_events = j_events_generator_.ComputeEvents(j_alignment);
        TRACE(ji);
        j_storages_.emplace_back(std::move(j_events));
        j_events_num += j_events.size();
    }
    TRACE(j_events_num << " events were computed for " << vdj_hits.JHitsNumber() << " J hits");
}

std::pair<recombination_utils::NongenomicInsertion, recombination_utils::NongenomicInsertion>
CustomHeavyChainRecombinationGenerator::RefineNongenomicInsertions(
        const recombination_utils::NongenomicInsertion &vd_insertion,
        const recombination_utils::NongenomicInsertion &dj_insertion) const
{
    recombination_utils::NongenomicInsertion vd_2(vd_insertion.StartPosition(), dj_insertion.EndPosition());
    recombination_utils::NongenomicInsertion dj_2(dj_insertion.EndPosition(), dj_insertion.EndPosition() - 1);
    return std::make_pair(vd_2, dj_2);
}

void CustomHeavyChainRecombinationGenerator::CreateRecombinations(
        recombination_utils::HcRecombinationStorage &recombination_storage,
        const recombination_utils::CleavedIgGeneAlignment &v_gene,
        const recombination_utils::CleavedIgGeneAlignment &d_gene,
        const recombination_utils::CleavedIgGeneAlignment &j_gene,
        recombination_utils::InsertionEventStorage &vd_insertions,
        recombination_utils::InsertionEventStorage &dj_insertions) const
{
    for(auto& vd_insertion : vd_insertions) {
        for(auto& dj_insertion: dj_insertions) {
            if(d_gene.GeneAlignment().Empty()) {
                auto refined_insertions = RefineNongenomicInsertions(vd_insertion, dj_insertion);
                vd_insertion = refined_insertions.first;
                dj_insertion = refined_insertions.second;
            }
            recombination_utils::HCRecombination recombination(recombination_storage.Read(),
                                                               v_gene, d_gene, j_gene,
                                                               vd_insertion, dj_insertion);
            TRACE(recombination);
            TRACE("-------");
            if(recombination.Valid())
                recombination_storage.AddRecombination(std::move(recombination));

        }
    }
    TRACE("#####");
}

void CustomHeavyChainRecombinationGenerator::CreateRecombinations(
        recombination_utils::HcRecombinationStorage &recombination_storage,
        const recombination_utils::IgGeneRecombinationEventStorage &v_events,
        const recombination_utils::IgGeneRecombinationEventStorage &d_events,
        const recombination_utils::IgGeneRecombinationEventStorage &j_events) const
{
    TRACE(v_events.size() << " V events were computed");
    TRACE(d_events.size() << " D events were computed");
    TRACE(j_events.size() << " J events were computed");
    for(auto vit = v_events.cbegin(); vit != v_events.cend(); vit++) {
        size_t v_end_position = (*vit).EndReadPosition();
        auto d_range = d_events.GetIndexRangeFromLeftPosition(v_end_position + 1);
        TRACE("D range: " << d_range.first << " - " << d_range.second);
        for(size_t di = d_range.first; di <= d_range.second; di++) {
            auto vd_insertions = vd_insertion_generator_.ComputeInsertionEvents(*vit, d_events[di]);
            size_t d_end_position = d_events[di].EndReadPosition();
            auto j_range = j_events.GetIndexRangeFromLeftPosition(d_end_position + 1);
            TRACE("J range: " << j_range.first << " - " << j_range.second);
            for(size_t ji = j_range.first; ji <= j_range.second; ji++) {
                auto dj_insertions = dj_insertion_generator_.ComputeInsertionEvents(d_events[di], j_events[ji]);
                CreateRecombinations(recombination_storage,
                                     *vit, d_events[di], j_events[ji],
                                     vd_insertions, dj_insertions);
            }
        }
    }
}

recombination_utils::HcRecombinationStorage CustomHeavyChainRecombinationGenerator::ComputeRecombinations(
    const VDJHits &vdj_hits)
{
    Clear();
    recombination_utils::HcRecombinationStorage recombination_storage(vdj_hits.ReadPtr());
    INFO("Generation of recombinations for read " << vdj_hits.Read().name);
    ComputeVEventStorages(vdj_hits);
    ComputeDEventStorages(vdj_hits);
    ComputeJEventStorages(vdj_hits);
    for(size_t vi = 0; vi < vdj_hits.VHitsNumber(); vi++) {
        // auto &v_alignment = vdj_hits.GetAlignmentByIndex(germline_utils::SegmentType::VariableSegment, vi);
        const recombination_utils::IgGeneRecombinationEventStorage &v_events = v_storages_[vi];
        for(size_t di = 0; di < vdj_hits.DHitsNumber(); di++) {
            // auto &d_alignment = vdj_hits.GetAlignmentByIndex(germline_utils::SegmentType::DiversitySegment, di);
            const recombination_utils::IgGeneRecombinationEventStorage &d_events = d_storages_[di];
            for(size_t ji = 0; ji < vdj_hits.JHitsNumber(); ji++) {
                // auto &j_alignment = vdj_hits.GetAlignmentByIndex(germline_utils::SegmentType::JoinSegment, ji);
                const recombination_utils::IgGeneRecombinationEventStorage &j_events = j_storages_[ji];
                // std::cout << "V. " << v_alignment.Alignment() << std::endl;
                // std::cout << "D. " << d_alignment.Alignment() << std::endl;
                // std::cout << "J. " << j_alignment.Alignment() << std::endl;
                CreateRecombinations(recombination_storage, v_events, d_events, j_events);
            }
        }
    }
    INFO(recombination_storage.size() << " valid recombinations were found");
    return recombination_storage;
}