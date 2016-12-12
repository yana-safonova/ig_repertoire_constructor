#pragma once

#include "../include_me.hpp"
#include "../gene_database.hpp"
#include "removing_settings.hpp"
#include "p_insertion_settings.hpp"
#include "n_insertion_settings.hpp"
#include "shm_settings.hpp"

// ----------------------------------------------------------------------------

class HC_VDJ_Recombination {
    const HC_GenesDatabase_Ptr db_ptr_;
    size_t vgene_index_;
    size_t dgene_index_;
    size_t jgene_index_;

    // endonuclease removing settings
    HC_RemovingSettings removing_settings_;

    // palindromic settings
    HC_PInsertionSettings p_insertion_settings_;

    // n-nucleotides insertion
    HC_NInsertionSettings n_insertion_settings_;

    // sequence and update flag
    string sequence_;
    bool update_sequence_;

    string ComputeSequence() const {
        return VgeneSeq().substr(0, VgeneLen()) +
                p_insertion_settings_.VEnd() + n_insertion_settings_.VD_Insertion() + p_insertion_settings_.DStart() +
                DgeneSeq().substr(removing_settings_.DStartLen(), DgeneLen()) +
                p_insertion_settings_.DEnd() + n_insertion_settings_.DJ_Insertion() + p_insertion_settings_.JStart() +
                JgeneSeq().substr(removing_settings_.JStartLen(), JgeneLen());
    }

public:
    HC_VDJ_Recombination(const HC_GenesDatabase_Ptr db_ptr,
            size_t vgene_index,
            size_t dgene_index,
            size_t jgene_index) :
            db_ptr_(db_ptr),
            vgene_index_(vgene_index),
            dgene_index_(dgene_index),
            jgene_index_(jgene_index),
            update_sequence_(true) {
    }

    size_t VgeneIndex() const { return vgene_index_; }

    size_t DgeneIndex() const { return dgene_index_; }

    size_t JgeneIndex() const { return jgene_index_; }

    // V gene block
    string VgeneName() const {
        return db_ptr_->GetByIndex(variable_gene, vgene_index_).GeneName();
    }

    string VshortGeneName() const {
        return db_ptr_->GetByIndex(variable_gene, vgene_index_).ShortGeneName();
    }

    string VgeneSeq() const {
        return db_ptr_->GetByIndex(variable_gene, vgene_index_).GeneSeq();
    }

    size_t VgeneLen() const {
        return VgeneSeq().size() - removing_settings_.VEndLen();
    }

    // D gene block
    string DgeneName() const {
        return db_ptr_->GetByIndex(diversity_gene, dgene_index_).GeneName();
    }

    string DshortGeneName() const {
        return db_ptr_->GetByIndex(diversity_gene, dgene_index_).ShortGeneName();
    }

    string DgeneSeq() const {
        return db_ptr_->GetByIndex(diversity_gene, dgene_index_).GeneSeq();
    }

    size_t DgeneLen() const {
        return DgeneSeq().size() - removing_settings_.DStartLen() - removing_settings_.DEndLen();
    }

    // D gene block
    string JgeneName() const {
        return db_ptr_->GetByIndex(join_gene, jgene_index_).GeneName();
    }

    string JgeneSeq() const {
        return db_ptr_->GetByIndex(join_gene, jgene_index_).GeneSeq();
    }

    string JshortGeneName() const {
        return db_ptr_->GetByIndex(join_gene, jgene_index_).ShortGeneName();
    }

    size_t JgeneLen() const {
        return JgeneSeq().size() - removing_settings_.JStartLen();
    }

    // settings get/set methods
    void AddRemovingSettings(HC_RemovingSettings removing_settings) {
        removing_settings_ = removing_settings;
        update_sequence_ = true;
    }

    const HC_RemovingSettings RemovingSettings() const { return removing_settings_; }

    void AddPInsertionSettings(HC_PInsertionSettings p_insertion_settings) {
        p_insertion_settings_ = p_insertion_settings;
        update_sequence_ = true;
    }

    const HC_PInsertionSettings PInsertionSettings() const { return p_insertion_settings_; }

    void AddNInsertionSettings(HC_NInsertionSettings n_insertion_settings) {
       n_insertion_settings_ = n_insertion_settings;
       update_sequence_ = true;
    }

    const HC_NInsertionSettings NInsertionSettings() const { return n_insertion_settings_; }

    const HC_GenesDatabase_Ptr GeneDB() const { return db_ptr_; }

    string Sequence() {
        if(update_sequence_)
            sequence_ = ComputeSequence();
        update_sequence_ = false;
        return sequence_;
    }

    typedef shared_ptr<HC_VDJ_Recombination> HC_VDJ_Recombination_Ptr;

    HC_VDJ_Recombination_Ptr Clone() {
        HC_VDJ_Recombination_Ptr vdj(new HC_VDJ_Recombination(db_ptr_, VgeneIndex(), DgeneIndex(), JgeneIndex()));
        vdj->AddRemovingSettings(RemovingSettings());
        vdj->AddNInsertionSettings(NInsertionSettings());
        vdj->AddPInsertionSettings(PInsertionSettings());
        return vdj;
    }

    string VJDRecombinationString() const {
        return VshortGeneName() + ";" + DshortGeneName() + ";" + JshortGeneName();
    }
};

ostream& operator<<(ostream &out, HC_VDJ_Recombination &obj) {
    out << "Indices: " << obj.VgeneIndex() << " " << obj.DgeneIndex() << " " << obj.JgeneIndex() << endl;
    out << "Vgene. " << obj.VgeneName() << endl;
    out << "Dgene. " << obj.DgeneName() << endl;
    out << "Jgene. " << obj.JgeneName() << endl;
    out << "Endonuclease removals:" << endl;
    out << obj.RemovingSettings();
    out << "P nucleotides settings:" << endl;
    out << obj.PInsertionSettings();
    out << "N nucleotides settings:" << endl;
    out << obj.NInsertionSettings();
    out << "Sequence: " << obj.Sequence() << endl;
    return out;
}

typedef HC_VDJ_Recombination::HC_VDJ_Recombination_Ptr HC_VDJ_Recombination_Ptr;

// ----------------------------------------------------------------------------

class LC_VDJ_Recombination {
    const LC_GenesDatabase_Ptr db_ptr_;
    size_t vgene_index_;
    size_t jgene_index_;

    // endonuclease removing settings
    LC_RemovingSettings removing_settings_;

    // palindromic settings
    LC_PInsertionSettings p_insertion_settings_;

    // n-nucleotides insertion
    LC_NInsertionSettings n_insertion_settings_;

    // sequence fields
    bool update_sequence_;
    string sequence_;

    string ComputeSequence() const {
        return VgeneSeq().substr(0, VgeneLen()) +
               p_insertion_settings_.VEnd() + n_insertion_settings_.VJ_Insertion() + p_insertion_settings_.JStart() +
               JgeneSeq().substr(removing_settings_.JStartLen(), JgeneLen());
    }

public:
    LC_VDJ_Recombination(const LC_GenesDatabase_Ptr db_ptr,
            size_t vgene_index,
            size_t jgene_index) :
            db_ptr_(db_ptr),
            vgene_index_(vgene_index),
            jgene_index_(jgene_index) { }

    // v block
    size_t VgeneIndex() const { return vgene_index_; }

    string VgeneSeq() const {
        return db_ptr_->GetByIndex(variable_gene, vgene_index_).GeneSeq();
    }

    string VgeneName() const {
        return db_ptr_->GetByIndex(variable_gene, vgene_index_).GeneName();
    }

    string VshortGeneName() const {
        return db_ptr_->GetByIndex(variable_gene, vgene_index_).ShortGeneName();
    }

    size_t VgeneLen() const {
        return VgeneSeq().size() - removing_settings_.VEndLen();
    }

    // j block
    size_t JgeneIndex() const { return jgene_index_; }

    string JgeneSeq() const {
        return db_ptr_->GetByIndex(join_gene, jgene_index_).GeneSeq();
    }

    string JgeneName() const {
        return db_ptr_->GetByIndex(join_gene, jgene_index_).GeneName();
    }

    string JshortGeneName() const {
        return db_ptr_->GetByIndex(join_gene, jgene_index_).ShortGeneName();
    }

    size_t JgeneLen() const {
        return JgeneSeq().size() - removing_settings_.JStartLen();
    }

    // removing settings
    void AddRemovingSettings(LC_RemovingSettings removing_settings) {
        removing_settings_ = removing_settings;
        update_sequence_ = true;
    }

    const LC_RemovingSettings RemovingSettings() const { return removing_settings_; }

    // p insertion settings
    void AddPInsertionSettings(LC_PInsertionSettings p_insertion_settings) {
        p_insertion_settings_ = p_insertion_settings;
        update_sequence_ = true;
    }

    const LC_PInsertionSettings PInsertionSettings() const { return p_insertion_settings_; }

    // n insertion settings
    void AddNInsertionSettings(LC_NInsertionSettings n_insertion_settings) {
        n_insertion_settings_ = n_insertion_settings;
        update_sequence_ = true;
    }

    const LC_NInsertionSettings NInsertionSettings() const { return n_insertion_settings_; }

    // db
    const LC_GenesDatabase_Ptr GeneDB() const { return db_ptr_; }

    string Sequence() {
        if(update_sequence_)
            sequence_ = ComputeSequence();
        update_sequence_ = false;
        return sequence_;
    }

    typedef shared_ptr<LC_VDJ_Recombination> LC_VDJ_Recombination_Ptr;

    LC_VDJ_Recombination_Ptr Clone() {
        LC_VDJ_Recombination_Ptr vdj(new LC_VDJ_Recombination(db_ptr_, VgeneIndex(), JgeneIndex()));
        vdj->AddRemovingSettings(RemovingSettings());
        vdj->AddNInsertionSettings(NInsertionSettings());
        vdj->AddPInsertionSettings(PInsertionSettings());
        return vdj;
    }

    string VJDRecombinationString() const {
        return VshortGeneName() + ";" + JshortGeneName();
    }
};

typedef LC_VDJ_Recombination::LC_VDJ_Recombination_Ptr LC_VDJ_Recombination_Ptr;

ostream& operator<<(ostream &out, const LC_VDJ_Recombination &obj) {
    out << "Indices: " << obj.VgeneIndex() << " " << obj.JgeneIndex() << endl;
    out << "Vgene. " << obj.VgeneName() << endl;
    out << "Jgene. " << obj.JgeneName() << endl;
    out << "Endonuclease removals:" << endl;
    out << obj.RemovingSettings();
    out << "P nucleotides settings:" << endl;
    out << obj.PInsertionSettings();
    out << "N nucleotides settings:" << endl;
    out << obj.NInsertionSettings();
    return out;
}

typedef shared_ptr<LC_VDJ_Recombination> LC_VDJ_Recombination_Ptr;