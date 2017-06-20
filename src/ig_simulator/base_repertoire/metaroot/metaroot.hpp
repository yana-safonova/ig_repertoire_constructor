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

    const annotation_utils::CDRLabeling cdr_labeling;

    std::string sequence;

    bool productive = true;

protected:
    static void PrepareGene(seqan::Dna5String& gene, int left_cleavage, int right_cleavage);
    virtual void CalculateSequence() = 0;
    virtual void print(std::ostream& out) const = 0;

public:
    AbstractMetaroot(const germline_utils::CustomGeneDatabase *v_db_p,
                     const germline_utils::CustomGeneDatabase *j_db_p,
                     const size_t v_ind,
                     const size_t j_ind,
                     const annotation_utils::CDRLabeling& cdr_labeling,
                     int cleavage_v,
                     int cleavage_j);

    const germline_utils::CustomGeneDatabase *V_DB_P() const { return v_db_p; }
    const germline_utils::CustomGeneDatabase *J_DB_P() const { return j_db_p; }

    size_t V_Ind()  const { return v_ind; }
    size_t J_Ind()  const { return j_ind; }
    int CleavageV() const { return cleavage_v; }
    int CleavageJ() const { return cleavage_j; }

    const annotation_utils::CDRLabeling CDRLabeling() const { return cdr_labeling; }

    size_t Length() const { return sequence.size(); }

    virtual const std::string& Sequence() const = 0;

    bool IsProductive() const { return productive; }
    void SetNonProductive()   { productive = false; }

    friend std::ostream& operator<<(std::ostream& out, const AbstractMetaroot& root) {
        root.print(out);
        return out;
    }

    AbstractMetaroot() = delete;
    AbstractMetaroot(const AbstractMetaroot&) = default;
    AbstractMetaroot(AbstractMetaroot&&) = default;
    AbstractMetaroot& operator=(const AbstractMetaroot&) = delete;
    AbstractMetaroot& operator=(AbstractMetaroot&&) = delete;

    virtual ~AbstractMetaroot() { }
};

using AbstractMetarootCPtr = std::unique_ptr<const AbstractMetaroot>;


class VJMetaroot final: public AbstractMetaroot {
private:
    const seqan::Dna5String insertion_vj;
    void CalculateSequence() override;
    void print(std::ostream&) const override;

public:
    VJMetaroot(const germline_utils::CustomGeneDatabase *v_db_p,
               const germline_utils::CustomGeneDatabase *j_db_p,
               const size_t v_ind,
               const size_t j_ind,
               const annotation_utils::CDRLabeling &cdr_labeling,
               int cleavage_v,
               int cleavage_j,
               seqan::Dna5String insertion_vj = "");

    const seqan::Dna5String& InsertionVJ() const { return insertion_vj; }

    const std::string& Sequence() const override;
};


class VDJMetaroot final: public AbstractMetaroot {
private:
    const germline_utils::CustomGeneDatabase * d_db_p;

    const size_t d_ind;

    // Negative cleavage means palindrome
    const int cleavage_d_left;
    const int cleavage_d_right;

    const seqan::Dna5String insertion_vd;
    const seqan::Dna5String insertion_dj;

private:
    void print(std::ostream&) const override;
    void CalculateSequence() override;

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
                const seqan::Dna5String& insertion_dj = "");

    const germline_utils::CustomGeneDatabase *D_DB_P() const { return d_db_p; }

    size_t D_Ind()                         const { return d_ind; }
    int CleavageDLeft()                    const { return cleavage_d_left; }
    int CleavageDRight()                   const { return cleavage_d_right; }
    const seqan::Dna5String& InsertionVD() const { return insertion_vd; }
    const seqan::Dna5String& InsertionDJ() const { return insertion_dj; }

    const std::string& Sequence() const override;
};

} // End namespace ig_simulator
