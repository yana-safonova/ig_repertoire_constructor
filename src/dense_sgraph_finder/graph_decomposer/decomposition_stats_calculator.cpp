#include "decomposition_stats_calculator.hpp"

using namespace dense_subgraph_finder;

size_t DecompositionStatsCalculator::GetClassOfVertices(size_t v1, size_t v2) {
    if (decomposition_->GetVertexClass(v1) == decomposition_->GetVertexClass(v2))
        return decomposition_->GetVertexClass(v1);
    return size_t(-1);
}

void DecompositionStatsCalculator::Initialize() {
    for (size_t i = 0; i < decomposition_->Size(); i++)
        num_edges_in_class.push_back(0);
    for (size_t i = 0; i < hamming_graph_->N(); i++)
        for (size_t j = hamming_graph_->RowIndex()[i]; j < hamming_graph_->RowIndex()[i + 1]; j++) {
            size_t v1 = i;
            size_t v2 = hamming_graph_->Col()[j];
            size_t class_id = GetClassOfVertices(v1, v2);
            if (class_id == size_t(-1))
                continue;
            num_edges_in_class[class_id]++;
        }
}

double DecompositionStatsCalculator::ComputeClassFillin(size_t class_id) {
    if (decomposition_->ClassSize(class_id) == 1)
        return 0;
    return double(num_edges_in_class[class_id]) / double(decomposition_->ClassSize(class_id) *
                                                         (decomposition_->ClassSize(class_id) - 1)) * 2;
}

size_t DecompositionStatsCalculator::ComputeClassSize(size_t class_id) {
    return decomposition_->ClassSize(class_id);
}

void DecompositionStatsCalculator::WriteStatsInFile(string filename) {
    TRACE("Initialization starts");
    Initialize();
    TRACE("Initialization ends");
    ofstream out(filename.c_str(), std::ios_base::app);
    for (size_t i = 0; i < decomposition_->Size(); i++) {
        assert(decomposition_->ClassSize(i) != 0);
        out << i << "\t" << decomposition_->ClassSize(i) << "\t" <<
        ComputeClassSize(i) << "\t" << ComputeClassFillin(i) << endl;
    }
    TRACE("Statistics were appended to " << filename);
    out.close();
}