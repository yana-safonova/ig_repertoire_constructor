#include "logger/logger.hpp"
#include "recombination_event_storage.hpp"

using namespace recombination_utils;

bool IgGeneRecombinationEventStorage::CheckConsistency(const CleavedIgGeneAlignment& new_gene_event) {
    if(recombination_events_.size() != 0) {
        size_t last_left_position = recombination_events_[recombination_events_.size() - 1].StartReadPosition();
        assert(last_left_position <= new_gene_event.StartReadPosition());
    }
    return new_gene_event.GeneAlignment()->subject().GeneType().Segment() == segment_type_;
}

void IgGeneRecombinationEventStorage::AddEvent(CleavedIgGeneAlignment new_gene_event) {
    if(CheckConsistency(new_gene_event))
        recombination_events_.push_back(new_gene_event);
    size_t left_position = new_gene_event.StartReadPosition();
    if(left_position_map_.find(left_position) == left_position_map_.end())
        left_position_map_[left_position] = Range(recombination_events_.size() - 1, size_t(-1));
    Range old_range = left_position_map_[left_position];
    left_position_map_[left_position] = Range(old_range.first, recombination_events_.size() - 1);
    min_left_position_ = std::min<size_t>(min_left_position_, new_gene_event.StartReadPosition());
    max_left_position_ = std::max<size_t>(max_left_position_, new_gene_event.StartReadPosition());
}

// both indices are inclusive
IgGeneRecombinationEventStorage::Range IgGeneRecombinationEventStorage::GetIndexRangeFromLeftPosition(size_t position) {
    if(position < min_left_position_)
        return Range(0, size() - 1);
    if(position > max_left_position_)
        return Range(size() + 1, size());
    return left_position_map_[position];
}

CleavedIgGeneAlignment IgGeneRecombinationEventStorage::operator[](size_t index) {
    assert(index < size());
    return recombination_events_[index];
}