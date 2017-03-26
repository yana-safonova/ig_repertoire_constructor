//
// Created by Andrew Bzikadze on 3/22/17.
//

#include "metaroot.hpp"
#include "convert.hpp"
#include <cstring>
#include "seqan/sequence.h"

namespace ig_simulator {

void AbstractMetaRoot::PrepareGene(seqan::Dna5String& gene, int left_cleavage, int right_cleavage) {
    VERIFY(static_cast<unsigned long>(abs(left_cleavage)) <= seqan::length(gene));
    if (left_cleavage > 0) {
        gene = seqan::suffix(gene, left_cleavage);
    } else if (left_cleavage < 0) {
        seqan::Dna5String pal = seqan::prefix(gene, -left_cleavage);
        seqan::reverseComplement(pal);
        pal += gene;
        gene = pal;
    }
    VERIFY(static_cast<unsigned long>(abs(right_cleavage)) <= seqan::length(gene));
    if (right_cleavage > 0) {
        gene = seqan::prefix(gene, seqan::length(gene) - right_cleavage);
    } else if (right_cleavage < 0) {
        seqan::Dna5String pal = seqan::suffix(gene, seqan::length(gene) + right_cleavage);
        seqan::reverseComplement(pal);
        gene += pal;
    }
}

const seqan::Dna5String& VJMetaRoot::Sequence() {
    if (sequence_calculated)
        return sequence;

    sequence_calculated = true;

    VERIFY(v_db_p != nullptr);
    VERIFY(j_db_p != nullptr);

    seqan::Dna5String v_gene = (*v_db_p)[v_ind].seq();
    seqan::Dna5String j_gene = (*j_db_p)[j_ind].seq();

    PrepareGene(v_gene, 0, cleavage_v);
    PrepareGene(j_gene, cleavage_j, 0);

    sequence = v_gene;
    sequence += insertion_vj;
    sequence += j_gene;

    return sequence;
}

const seqan::Dna5String& VDJMetaRoot::Sequence() {
    if (sequence_calculated)
        return sequence;

    sequence_calculated = true;

    VERIFY(v_db_p != nullptr);
    VERIFY(d_db_p != nullptr);
    VERIFY(j_db_p != nullptr);

    seqan::Dna5String v_gene = (*v_db_p)[v_ind].seq();
    seqan::Dna5String d_gene = (*d_db_p)[d_ind].seq();
    seqan::Dna5String j_gene = (*j_db_p)[j_ind].seq();

    PrepareGene(v_gene, 0, cleavage_v);
    PrepareGene(d_gene, cleavage_d_left, cleavage_d_right);
    PrepareGene(j_gene, cleavage_j, 0);

    sequence = v_gene;
    sequence += insertion_vd;
    sequence += d_gene;
    sequence += insertion_dj;
    sequence += j_gene;

    return sequence;
}

std::ostream& operator<<(std::ostream& out, const VJMetaRoot& root) {
    out << "VJMetaroot:\n\n" <<

        "Index (from 1) of V in database " << root.V_Ind() + 1 << " / " << root.V_DB_P()->size() << "\n" <<
        "Index (from 1) of J in database " << root.J_Ind() + 1 << " / " << root.J_DB_P()->size() << "\n\n" <<

        "V gene: " << (*root.V_DB_P())[root.V_Ind()].seq() << "\n" <<
        "J gene: " << (*root.J_DB_P())[root.J_Ind()].seq() << "\n\n" <<

        "Cleavage in V gene: " << root.CleavageV() << "\n" <<
        "Cleavage in J gene: " << root.CleavageJ() << "\n\n" <<

        "Insertion in VJ junction: " << root.InsertionVJ() << "\n\n" <<

        "CDR1: " << root.CDRLabeling().cdr1.start_pos << " " << root.CDRLabeling().cdr1.end_pos << "\n" <<
        "CDR2: " << root.CDRLabeling().cdr2.start_pos << " " << root.CDRLabeling().cdr2.end_pos << "\n" <<
        "CDR3: " << root.CDRLabeling().cdr3.start_pos << " " << root.CDRLabeling().cdr3.end_pos << "\n";

    return out;
}

std::ostream& operator<<(std::ostream& out, const VDJMetaRoot& root) {
    out << "VDJMetaroot:\n\n" <<

        "Index (from 1) of V in database " << root.V_Ind() + 1 << " / " << root.V_DB_P()->size() << "\n" <<
        "Index (from 1) of D in database " << root.D_Ind() + 1 << " / " << root.D_DB_P()->size() << "\n" <<
        "Index (from 1) of J in database " << root.J_Ind() + 1 << " / " << root.J_DB_P()->size() << "\n\n" <<

        "V gene: " << (*root.V_DB_P())[root.V_Ind()].seq() << "\n" <<
        "D gene: " << (*root.D_DB_P())[root.D_Ind()].seq() << "\n" <<
        "J gene: " << (*root.J_DB_P())[root.J_Ind()].seq() << "\n\n" <<

        "Cleavage in V gene: " << root.CleavageV() << "\n" <<
        "Cleavage in D gene (left): " << root.CleavageDLeft() << "\n" <<
        "Cleavage in D gene (right): " << root.CleavageDRight() << "\n" <<
        "Cleavage in J gene: " << root.CleavageJ() << "\n\n" <<

        "Insertion in VD junction: " << root.InsertionVD() << "\n" <<
        "Insertion in DJ junction: " << root.InsertionDJ() << "\n\n" <<

        "CDR1: " << root.CDRLabeling().cdr1.start_pos << " " << root.CDRLabeling().cdr1.end_pos << "\n" <<
        "CDR2: " << root.CDRLabeling().cdr2.start_pos << " " << root.CDRLabeling().cdr2.end_pos << "\n" <<
        "CDR3: " << root.CDRLabeling().cdr3.start_pos << " " << root.CDRLabeling().cdr3.end_pos << "\n";
    return out;
}

} // End namespace ig_simulator
