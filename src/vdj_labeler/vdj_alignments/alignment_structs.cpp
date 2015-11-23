#include "alignment_structs.hpp"

#include "seqan/sequence.h"
#include <seqan/stream.h>

using namespace std;

ostream& operator<<(ostream &out, const AlignmentPositions& obj) {
    out << "Query pos: " << obj.query_pos.first << " - " << obj.query_pos.second << ". Subject pos: " <<
            obj.subject_pos.first << " - " << obj.subject_pos.second;
    return out;
}

ostream& operator<<(ostream& out, const IgGeneAlignmentPositions& obj) {
    out << "Read. " << *obj.read << endl;
    out << "Gene. " << *obj.ig_gene << endl;
    out << obj.alignment << endl;
    return out;
}

ostream& operator<<(ostream &out, const IgGeneAlignment& ig_gene_alignment) {
    out << ig_gene_alignment.positions << endl;
    out << ig_gene_alignment.alignment;
    return out;
}