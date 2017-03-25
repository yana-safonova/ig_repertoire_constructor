//
// Created by Andrew Bzikadze on 3/22/17.
//

#pragma once

#include "gtest/gtest_prod.h"
#include "germline_utils/germline_db_generator.hpp"

namespace ig_simulator {

class AbstractMetaRoot {
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

    bool sequence_calculated;
    seqan::Dna5String sequence;

    static void PrepareGene(seqan::Dna5String& gene, int left_cleavage, int right_cleavage);

public:
    AbstractMetaRoot(const germline_utils::CustomGeneDatabase *v_db_p = nullptr,
                     const germline_utils::CustomGeneDatabase *j_db_p = nullptr,
                     const size_t v_ind = size_t(-1),
                     const size_t j_ind = size_t(-1),
                     int cleavage_v = 0,
                     int cleavage_j = 0):
            v_db_p(v_db_p),
            j_db_p(j_db_p),
            v_ind(v_ind),
            j_ind(j_ind),
            cleavage_v(cleavage_v),
            cleavage_j(cleavage_j),
            sequence_calculated(false) {
        if (v_db_p != nullptr) {
            VERIFY(v_ind < v_db_p->size());
        } else {
            VERIFY(v_ind == size_t(-1));
            VERIFY(cleavage_v == 0);
        }

        if (j_db_p != nullptr) {
            VERIFY(j_ind < j_db_p->size());
        } else {
            VERIFY(j_ind == size_t(-1));
            VERIFY(cleavage_j == 0);
        }
    }

    const germline_utils::CustomGeneDatabase *V_DB_P() const { return v_db_p; }
    const germline_utils::CustomGeneDatabase *J_DB_P() const { return j_db_p; }

    size_t V_Ind()  const { return v_ind; }
    size_t J_Ind()  const { return j_ind; }
    int CleavageV() const { return cleavage_v; }
    int CleavageJ() const { return cleavage_j; }

    virtual const seqan::Dna5String& Sequence() = 0;
    virtual ~AbstractMetaRoot() { }
};


class VJMetaRoot : public AbstractMetaRoot {
private:
    const seqan::Dna5String insertion_vj;

public:
    VJMetaRoot(const germline_utils::CustomGeneDatabase *v_db_p = nullptr,
               const germline_utils::CustomGeneDatabase *j_db_p = nullptr,
               const size_t v_ind = size_t(-1),
               const size_t j_ind = size_t(-1),
               int cleavage_v = 0,
               int cleavage_j = 0,
               seqan::Dna5String insertion_vj = "") :
        AbstractMetaRoot(v_db_p, j_db_p, v_ind, j_ind, cleavage_v, cleavage_j),
        insertion_vj(insertion_vj)
    { }

    const seqan::Dna5String& InsertionVJ() const { return insertion_vj; }

    virtual const seqan::Dna5String& Sequence() override;
};

class VDJMetaRoot : public AbstractMetaRoot {
private:
    const germline_utils::CustomGeneDatabase * d_db_p;

    const size_t d_ind;

    // Negative cleavage means palindrome
    const int cleavage_d_left;
    const int cleavage_d_right;

    const seqan::Dna5String insertion_vd;
    const seqan::Dna5String insertion_dj;

public:
    VDJMetaRoot(const germline_utils::CustomGeneDatabase *v_db_p = nullptr,
                const germline_utils::CustomGeneDatabase *d_db_p = nullptr,
                const germline_utils::CustomGeneDatabase *j_db_p = nullptr,
                const size_t v_ind = size_t(-1),
                const size_t d_ind = size_t(-1),
                const size_t j_ind = size_t(-1),
                int cleavage_v = 0,
                int cleavage_d_left = 0,
                int cleavage_d_right = 0,
                int cleavage_j = 0,
                const seqan::Dna5String& insertion_vd = "",
                const seqan::Dna5String& insertion_dj = "") :
            AbstractMetaRoot(v_db_p, j_db_p, v_ind, j_ind, cleavage_v, cleavage_j),
            d_db_p(d_db_p),
            d_ind(d_ind),
            cleavage_d_left(cleavage_d_left),
            cleavage_d_right(cleavage_d_right),
            insertion_vd(insertion_vd),
            insertion_dj(insertion_dj) {
        if (d_db_p != nullptr) {
            VERIFY(d_ind < d_db_p->size());
        } else {
            VERIFY(d_ind == size_t(-1));
            VERIFY(cleavage_j == 0);
        }
    }

    const germline_utils::CustomGeneDatabase *D_DB_P() const { return d_db_p; }

    size_t D_Ind()                         const { return d_ind; }
    int CleavageDLeft()                    const { return cleavage_d_left; }
    int CleavageDRight()                   const { return cleavage_d_right; }
    const seqan::Dna5String& InsertionVD() const { return insertion_vd; }
    const seqan::Dna5String& InsertionDJ() const { return insertion_dj; }

    virtual const seqan::Dna5String& Sequence() override;
};

std::ostream& operator<<(std::ostream& out, const VJMetaRoot& root);
std::ostream& operator<<(std::ostream& out, const VDJMetaRoot& root);

} // End namespace ig_simulator
