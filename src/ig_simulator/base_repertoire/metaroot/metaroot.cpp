//
// Created by Andrew Bzikadze on 3/22/17.
//

#include "metaroot.hpp"
#include "convert.hpp"
#include <cstring>
#include "seqan/sequence.h"

namespace ig_simulator {

inline
AbstractMetaroot::AbstractMetaroot(const germline_utils::CustomGeneDatabase *v_db_p,
                                   const germline_utils::CustomGeneDatabase *j_db_p,
                                   const size_t v_ind,
                                   const size_t j_ind,
                                   const annotation_utils::CDRLabeling& cdr_labeling,
                                   int cleavage_v,
                                   int cleavage_j) :
    v_db_p(check_pointer(v_db_p)),
    j_db_p(check_pointer(j_db_p)),
    v_ind(v_ind),
    j_ind(j_ind),
    cleavage_v(cleavage_v),
    cleavage_j(cleavage_j),
    cdr_labeling(cdr_labeling)
{
    VERIFY(v_ind < v_db_p->size());
    VERIFY(j_ind < j_db_p->size());
}

void AbstractMetaroot::PrepareGene(seqan::Dna5String& gene, int left_cleavage, int right_cleavage) {
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

VJMetaroot::VJMetaroot(const germline_utils::CustomGeneDatabase *v_db_p,
                       const germline_utils::CustomGeneDatabase *j_db_p,
                       const size_t v_ind,
                       const size_t j_ind,
                       const annotation_utils::CDRLabeling &cdr_labeling,
                       int cleavage_v,
                       int cleavage_j,
                       seqan::Dna5String insertion_vj) :
    AbstractMetaroot(v_db_p, j_db_p, v_ind, j_ind, cdr_labeling, cleavage_v, cleavage_j),
    insertion_vj(insertion_vj)
{
    CalculateSequence();
}

const std::string& VJMetaroot::Sequence() const { return sequence; }

void VJMetaroot::CalculateSequence() {
    VERIFY(v_db_p != nullptr);
    VERIFY(j_db_p != nullptr);

    seqan::Dna5String v_gene = (*v_db_p)[v_ind].seq();
    seqan::Dna5String j_gene = (*j_db_p)[j_ind].seq();

    PrepareGene(v_gene, 0, cleavage_v);
    PrepareGene(j_gene, cleavage_j, 0);

    sequence = core::seqan_string_to_string(v_gene);
    sequence += core::seqan_string_to_string(insertion_vj);
    sequence += core::seqan_string_to_string(j_gene);
}

VDJMetaroot::VDJMetaroot(const germline_utils::CustomGeneDatabase *v_db_p,
                         const germline_utils::CustomGeneDatabase *d_db_p,
                         const germline_utils::CustomGeneDatabase *j_db_p,
                         const size_t v_ind,
                         const size_t d_ind,
                         const size_t j_ind,
                         const annotation_utils::CDRLabeling& cdr_labeling,
                         int cleavage_v,
                         int cleavage_d_left,
                         int cleavage_d_right,
                         int cleavage_j,
                         const seqan::Dna5String& insertion_vd,
                         const seqan::Dna5String& insertion_dj) :
    AbstractMetaroot(v_db_p, j_db_p, v_ind, j_ind, cdr_labeling, cleavage_v, cleavage_j),
    d_db_p(check_pointer(d_db_p)),
    d_ind(d_ind),
    cleavage_d_left(cleavage_d_left),
    cleavage_d_right(cleavage_d_right),
    insertion_vd(insertion_vd),
    insertion_dj(insertion_dj)
{
    VERIFY(d_ind < d_db_p->size());
    CalculateSequence();
}

const std::string& VDJMetaroot::Sequence() const { return sequence; }

void VDJMetaroot::CalculateSequence() {
    VERIFY(v_db_p != nullptr);
    VERIFY(d_db_p != nullptr);
    VERIFY(j_db_p != nullptr);

    seqan::Dna5String v_gene = (*v_db_p)[v_ind].seq();
    seqan::Dna5String d_gene = (*d_db_p)[d_ind].seq();
    seqan::Dna5String j_gene = (*j_db_p)[j_ind].seq();

    PrepareGene(v_gene, 0, cleavage_v);
    PrepareGene(d_gene, cleavage_d_left, cleavage_d_right);
    PrepareGene(j_gene, cleavage_j, 0);

    sequence = core::seqan_string_to_string(v_gene);
    sequence += core::seqan_string_to_string(insertion_vd);
    sequence += core::seqan_string_to_string(d_gene);
    sequence += core::seqan_string_to_string(insertion_dj);
    sequence += core::seqan_string_to_string(j_gene);
}

void VJMetaroot::print(std::ostream& out) const {
    out << "VJMetaroot:\n\n" <<

        "V gene name: " << (*V_DB_P())[V_Ind()].name() << "\n" <<
        "J gene name: " << (*J_DB_P())[J_Ind()].name() << "\n\n" <<

        "V gene: " << (*V_DB_P())[V_Ind()].seq() << "\n" <<
        "J gene: " << (*J_DB_P())[J_Ind()].seq() << "\n\n" <<

        "Cleavage in V gene: " << CleavageV() << "\n" <<
        "Cleavage in J gene: " << CleavageJ() << "\n\n" <<

        "Insertion in VJ junction: " << InsertionVJ() << "\n\n" <<

        "CDR1 positions: " << CDRLabeling().cdr1.start_pos << " " << CDRLabeling().cdr1.end_pos << "\n" <<
        "CDR2 positions: " << CDRLabeling().cdr2.start_pos << " " << CDRLabeling().cdr2.end_pos << "\n" <<
        "CDR3 positions: " << CDRLabeling().cdr3.start_pos << " " << CDRLabeling().cdr3.end_pos << "\n\n" <<

        "CDR1: " << sequence.substr(CDRLabeling().cdr1.start_pos,
                                    CDRLabeling().cdr1.end_pos - CDRLabeling().cdr1.start_pos + 1) << "\n" <<
        "CDR2: " << sequence.substr(CDRLabeling().cdr2.start_pos,
                                    CDRLabeling().cdr2.end_pos - CDRLabeling().cdr2.start_pos + 1) << "\n" <<
        "CDR3: " << sequence.substr(CDRLabeling().cdr3.start_pos,
                                    CDRLabeling().cdr3.end_pos - CDRLabeling().cdr3.start_pos + 1) << "\n";
}

void VDJMetaroot::print(std::ostream& out) const {
    out << "VDJMetaroot:\n\n" <<

        "V gene name: " << (*V_DB_P())[V_Ind()].name() << "\n" <<
        "D gene name: " << (*D_DB_P())[D_Ind()].name() << "\n" <<
        "J gene name: " << (*J_DB_P())[J_Ind()].name() << "\n\n" <<

        "V gene: " << (*V_DB_P())[V_Ind()].seq() << "\n" <<
        "D gene: " << (*D_DB_P())[D_Ind()].seq() << "\n" <<
        "J gene: " << (*J_DB_P())[J_Ind()].seq() << "\n\n" <<

        "Cleavage in V gene: " << CleavageV() << "\n" <<
        "Cleavage in D gene (left): " << CleavageDLeft() << "\n" <<
        "Cleavage in D gene (right): " << CleavageDRight() << "\n" <<
        "Cleavage in J gene: " << CleavageJ() << "\n\n" <<

        "Insertion in VD junction: " << InsertionVD() << "\n" <<
        "Insertion in DJ junction: " << InsertionDJ() << "\n\n" <<

        "CDR1 positions: " << CDRLabeling().cdr1.start_pos << " " << CDRLabeling().cdr1.end_pos << "\n" <<
        "CDR2 positions: " << CDRLabeling().cdr2.start_pos << " " << CDRLabeling().cdr2.end_pos << "\n" <<
        "CDR3 positions: " << CDRLabeling().cdr3.start_pos << " " << CDRLabeling().cdr3.end_pos << "\n\n" <<


        "CDR1: " << sequence.substr(CDRLabeling().cdr1.start_pos,
                                    CDRLabeling().cdr1.end_pos - CDRLabeling().cdr1.start_pos + 1) << "\n" <<
        "CDR2: " << sequence.substr(CDRLabeling().cdr2.start_pos,
                                    CDRLabeling().cdr2.end_pos - CDRLabeling().cdr2.start_pos + 1) << "\n" <<
        "CDR3: " << sequence.substr(CDRLabeling().cdr3.start_pos,
                                    CDRLabeling().cdr3.end_pos - CDRLabeling().cdr3.start_pos + 1) << "\n";
}

} // End namespace ig_simulator
