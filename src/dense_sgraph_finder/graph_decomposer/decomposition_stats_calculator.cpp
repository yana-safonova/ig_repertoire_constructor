#include "decomposition_stats_calculator.hpp"

using namespace dense_subgraph_finder;

size_t DecompositionStatsCalculator::GetClassOfVertices(size_t v1, size_t v2) {
    if (decomposition_->GetVertexClass(v1) == decomposition_->GetVertexClass(v2))
        return decomposition_->GetVertexClass(v1);
    return size_t(-1);
}

void DecompositionStatsCalculator::Initialize() {
    for (size_t i = 0; i < decomposition_->Size(); i++) {
        num_edges_in_class_.push_back(0);
        class_edge_fillin_.push_back(0.0);
    }
    for (size_t i = 0; i < hamming_graph_->N(); i++)
        for (size_t j = hamming_graph_->RowIndex()[i]; j < hamming_graph_->RowIndex()[i + 1]; j++) {
            size_t v1 = i;
            size_t v2 = hamming_graph_->Col()[j];
            size_t class_id = GetClassOfVertices(v1, v2);
            if (class_id == size_t(-1))
                continue;
            num_edges_in_class_[class_id]++;
        }
    for (size_t i = 0; i < decomposition_->Size(); i++) {
        class_edge_fillin_[i] = ComputeClassFillin(i);
    }
    ComputeShortStats();
}

void DecompositionStatsCalculator::ComputeShortStats() {
    size_t max_class_size = 0;
    size_t max_class_id = size_t(-1);
    size_t num_nontrivial_classes = 0;
    for(size_t i = 0; i < class_edge_fillin_.size(); i++) {
        if(decomposition_->ClassSize(i) > 1) {
            max_edge_fillin_ = max<double>(class_edge_fillin_[i], max_edge_fillin_);
            average_fillin_ += class_edge_fillin_[i];
            num_nontrivial_classes++;
        }
        else num_trivial_classes_++;
        if(decomposition_->ClassSize(i) > max_class_size) {
            max_class_size = decomposition_->ClassSize(i);
            max_class_id = i;
        }
    }
    average_fillin_ /= double(num_nontrivial_classes);
    if(num_nontrivial_classes == 0)
        average_fillin_ = 0.0;
    if(max_class_id == size_t(-1))
        fillin_of_max_class_ = 0.0;
    else
        fillin_of_max_class_ = class_edge_fillin_[max_class_id];
}

double DecompositionStatsCalculator::ComputeClassFillin(size_t class_id) {
    if (decomposition_->ClassSize(class_id) == 1)
        return 0;
    return double(num_edges_in_class_[class_id]) / double(decomposition_->ClassSize(class_id) *
                                                         (decomposition_->ClassSize(class_id) - 1)) * 2;
}

void DecompositionStatsCalculator::WriteAllStats(ostream &out) {
    TRACE("Initialization starts");
    Initialize();
    TRACE("Initialization ends");
    for (size_t i = 0; i < decomposition_->Size(); i++) {
        assert(decomposition_->ClassSize(i) != 0);
        out << i << "\t" << class_edge_fillin_[i] << endl;
    }
}

void DecompositionStatsCalculator::WriteShortStats(ostream &) {
    INFO("Decomposition contains " << decomposition_->Size() << " classes");
    float trivial_classes_perc = float(num_trivial_classes_ * 100) / float(decomposition_->Size());
    if(decomposition_->Size() == 0)
        trivial_classes_perc = 0.0;
    INFO("# trivial classes: " << num_trivial_classes_ << " (" << trivial_classes_perc << "%)");
    INFO("Size of maximal class: " << decomposition_->MaxClassSize());
    INFO("Edge fillin of maximal class: " << fillin_of_max_class_);
    INFO("Maximal edge fillin: " << max_edge_fillin_);
    INFO("Average edge fillin: " << average_fillin_);
}