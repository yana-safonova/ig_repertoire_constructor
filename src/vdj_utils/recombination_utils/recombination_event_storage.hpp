#pragma once

#include "cleaved_gene.hpp"
#include <unordered_map>
#include "germline_utils/germline_gene_type.hpp"
#include "cleaved_gene.hpp"

namespace recombination_utils {

class IgGeneRecombinationEventStorage {
    germline_utils::SegmentType segment_type_;
    std::vector <CleavedIgGeneAlignment> recombination_events_;
    std::unordered_map <size_t, std::pair<size_t, size_t>> left_position_map_;
    size_t min_left_position_;
    size_t max_left_position_;

    // storage consider that events are added in the order of left event length decreasing:
    // from the maximal cleavage to maximal palindrome
    bool CheckConsistency(const CleavedIgGeneAlignment &new_gene_event);

public:
    IgGeneRecombinationEventStorage(germline_utils::SegmentType gene_type) :
        segment_type_(gene_type) {
        min_left_position_ = size_t(-1);
        max_left_position_ = 0;
    }

    void AddEvent(CleavedIgGeneAlignment new_gene_event);

    typedef std::vector<CleavedIgGeneAlignment>::const_iterator gene_segment_event_iterator;

    gene_segment_event_iterator cbegin() const { return recombination_events_.cbegin(); }

    gene_segment_event_iterator cend() const { return recombination_events_.cend(); }

    size_t size() const { return recombination_events_.size(); }

    typedef std::pair <size_t, size_t> Range;

    Range GetIndexRangeFromLeftPosition(size_t position);

    CleavedIgGeneAlignment operator[](size_t index);
};

typedef std::shared_ptr <IgGeneRecombinationEventStorage> IgGeneRecombinationEventStoragePtr;

} // End namespace vdj_labeler