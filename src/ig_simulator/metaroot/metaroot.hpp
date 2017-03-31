//
// Created by Andrew Bzikadze on 3/22/17.
//

#pragma once

#include "gtest/gtest_prod.h"
#include "germline_utils/germline_db_generator.hpp"
#include "annotation_utils/cdr_labeling_primitives.hpp"
#include "ig_simulator_utils.hpp"

namespace ig_simulator {

class AbstractMetaroot {
    friend class IgSimulatorTest;
    FRIEND_TEST(IgSimulatorTest, PrepareGeneTest);

protected:
    const germline_utils::CustomGeneDatabase * v_db_p;
    const germline_utils::CustomGeneDatabase * j_db_p;

    const size_t v_ind;
    const size_t j_ind;

    // Negative cleavage means palindrome
    const int cleavage_v;
    const int cleavage_j;

    annotation_utils::CDRLabeling cdr_labeling;

    seqan::Dna5String sequence;

    static void PrepareGene(seqan::Dna5String& gene, int left_cleavage, int right_cleavage);
    virtual void CalculateSequence() = 0;

public:
    AbstractMetaroot(const germline_utils::CustomGeneDatabase *v_db_p,
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

    const germline_utils::CustomGeneDatabase *V_DB_P() const { return v_db_p; }
    const germline_utils::CustomGeneDatabase *J_DB_P() const { return j_db_p; }

    size_t V_Ind()  const { return v_ind; }
    size_t J_Ind()  const { return j_ind; }
    int CleavageV() const { return cleavage_v; }
    int CleavageJ() const { return cleavage_j; }

    const annotation_utils::CDRLabeling CDRLabeling() const { return cdr_labeling; }

    virtual const seqan::Dna5String& Sequence() const = 0;
    virtual ~AbstractMetaroot() { }
};

using AbstractMetarootPtr = std::unique_ptr<AbstractMetaroot>;


class VJMetaroot : public AbstractMetaroot {
private:
    const seqan::Dna5String insertion_vj;
    virtual void CalculateSequence() override;

public:
    VJMetaroot(const germline_utils::CustomGeneDatabase *v_db_p,
               const germline_utils::CustomGeneDatabase *j_db_p,
               const size_t v_ind,
               const size_t j_ind,
               const annotation_utils::CDRLabeling &cdr_labeling,
               int cleavage_v,
               int cleavage_j,
               seqan::Dna5String insertion_vj = "") :
        AbstractMetaroot(v_db_p, j_db_p, v_ind, j_ind, cdr_labeling, cleavage_v, cleavage_j),
        insertion_vj(insertion_vj)
    {
        CalculateSequence();
    }

    const seqan::Dna5String& InsertionVJ() const { return insertion_vj; }

    virtual const seqan::Dna5String& Sequence() const override;
};

class VDJMetaroot : public AbstractMetaroot {
private:
    const germline_utils::CustomGeneDatabase * d_db_p;

    const size_t d_ind;

    // Negative cleavage means palindrome
    const int cleavage_d_left;
    const int cleavage_d_right;

    const seqan::Dna5String insertion_vd;
    const seqan::Dna5String insertion_dj;

private:
    virtual void CalculateSequence() override;

public:
    VDJMetaroot(const germline_utils::CustomGeneDatabase *v_db_p,
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
                const seqan::Dna5String& insertion_vd = "",
                const seqan::Dna5String& insertion_dj = "") :
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

    const germline_utils::CustomGeneDatabase *D_DB_P() const { return d_db_p; }

    size_t D_Ind()                         const { return d_ind; }
    int CleavageDLeft()                    const { return cleavage_d_left; }
    int CleavageDRight()                   const { return cleavage_d_right; }
    const seqan::Dna5String& InsertionVD() const { return insertion_vd; }
    const seqan::Dna5String& InsertionDJ() const { return insertion_dj; }

    virtual const seqan::Dna5String& Sequence() const override;
};

std::ostream& operator<<(std::ostream& out, const VJMetaroot& root);
std::ostream& operator<<(std::ostream& out, const VDJMetaroot& root);

} // End namespace ig_simulator
