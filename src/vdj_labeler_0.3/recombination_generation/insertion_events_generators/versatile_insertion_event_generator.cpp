#include "versatile_insertion_event_generator.hpp"

using namespace vdj_labeler;

recombination_utils::InsertionEventStoragePtr VersatileInsertionGenerator::ComputeInsertionEvents(
        recombination_utils::CleavedIgGeneAlignment left_gene_alignment,
        recombination_utils::CleavedIgGeneAlignment right_gene_alignment)
{
    recombination_utils::InsertionEventStoragePtr insertions(new recombination_utils::InsertionEventStorage());
    recombination_utils::NongenomicInsertion insertion(
        left_gene_alignment.EndReadPosition() + 1,
        right_gene_alignment.StartReadPosition() - 1);
    insertions->AddInsertion(insertion);
    return insertions;
}