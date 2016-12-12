#pragma once

#include "ig_structs/ig_structure_structs.hpp"
#include "ig_structs/shm_settings.hpp"
#include "ig_structs/cdr_settings.hpp"

template<class VDJ_Recombination_Ptr>
class IgVariableRegion {
    VDJ_Recombination_Ptr vdj_recombination_;
    SHMSettings shm_settings_;
    CDRSettings cdr_settings_;

    string sequence_;
    bool update_seq_;

    void ComputeSequence() {
        string new_sequence = "";
        size_t old_seq_ptr = 0;
        for(auto it = shm_settings_.begin(); it != shm_settings_.end(); it++) {
            size_t shm_pos = it->first;
            // append all between old ptr and shm pos
            for(size_t i = old_seq_ptr; i < shm_pos; i++)
                new_sequence += sequence_[i];
            if(it->second.IsSubstitution()) {
                new_sequence += it->second.GetSubstitution();
                old_seq_ptr = shm_pos + 1;
            }
            else if(it->second.IsDeletion()) {
                old_seq_ptr = shm_pos + it->second.GetNumDeleted();
            }
            else if(it->second.IsInsertion()) {
                new_sequence += it->second.GetInsertedString();
                old_seq_ptr = shm_pos + 1;
            }
        }
        for(size_t i = old_seq_ptr; i < sequence_.size(); i++)
            new_sequence += sequence_[i];
        sequence_ = new_sequence;
        update_seq_ = false;
    }

public:
    IgVariableRegion(VDJ_Recombination_Ptr vdj_recombination) :
            vdj_recombination_(vdj_recombination),
            sequence_(vdj_recombination->Sequence()),
            update_seq_(false) {
    }

    VDJ_Recombination_Ptr VDJ_Recombination() {
        return vdj_recombination_;
    }

    void SetSHMSettings(SHMSettings shm_settings) {
        shm_settings_ = shm_settings;
        update_seq_ = true;
    }

    const SHMSettings& GetSHMSettings() const { return shm_settings_; };

    void SetCDRSettings(CDRSettings cdr_settings) {
        cdr_settings_ = cdr_settings;
    }

    const CDRSettings& GetCDRSettings() const { return cdr_settings_; }

    string Sequence() {
        if(update_seq_)
            ComputeSequence();
        return sequence_;
    }

    size_t Length() {
        return Sequence().size();
    }

    typedef shared_ptr<IgVariableRegion<VDJ_Recombination_Ptr> > IgVariableRegionPtr;

    IgVariableRegionPtr Clone() {
        IgVariableRegionPtr variable_region_ptr(
                new IgVariableRegion<VDJ_Recombination_Ptr>(vdj_recombination_->Clone()));
        variable_region_ptr->SetSHMSettings(shm_settings_);
        variable_region_ptr->SetCDRSettings(cdr_settings_);
        return variable_region_ptr;
    }
};

typedef IgVariableRegion<HC_VDJ_Recombination_Ptr> HC_VariableRegion;
typedef IgVariableRegion<LC_VDJ_Recombination_Ptr> LC_VariableRegion;
typedef shared_ptr<IgVariableRegion<HC_VDJ_Recombination_Ptr> > HC_VariableRegionPtr;
typedef shared_ptr<IgVariableRegion<LC_VDJ_Recombination_Ptr> > LC_VariableRegionPtr;

// ----------------------------------------------------------------------------

// cluster is characterized by:
//      - IgVariableRegion (VDJ recombination, CDR3 labeling)
//      - multiplicity

template<class IgVariableRegionPtr>
class IgCluster {
    IgVariableRegionPtr variable_region_;
    size_t multiplicity_;

public:
    IgCluster(IgVariableRegionPtr variable_region, size_t multiplicity) :
            variable_region_(variable_region),
            multiplicity_(multiplicity) { }

    IgVariableRegionPtr IgVariableRegion() const { return variable_region_; }

    size_t Multiplicity() const { return multiplicity_; }

    string Sequence() const { return variable_region_->Sequence(); }
};

typedef IgCluster<HC_VariableRegionPtr> HC_Cluster;
typedef IgCluster<LC_VariableRegionPtr> LC_Cluster;

typedef vector<HC_Cluster>::const_iterator HC_ClusterIterator;
typedef vector<LC_Cluster>::const_iterator LC_ClusterIterator;

// ----------------------------------------------------------------------------

template<class IgCluster, class IgClusterIterator>
class Repertoire {
    vector<IgCluster> ig_clusters_;

public:
    Repertoire() { }

    void Add(IgCluster ig_cluster, size_t cluster_multiplicity = 1) {
        for(size_t i = 0; i < cluster_multiplicity; i++)
            ig_clusters_.push_back(ig_cluster);
    }

    size_t Size() const {
        return ig_clusters_.size();
    }

    IgClusterIterator begin() const {
        return ig_clusters_.begin();
    }

    IgClusterIterator end() const {
        return ig_clusters_.end();
    }

    size_t NumberAntibodies() {
        size_t size = 0;
        for(auto it = begin(); it != end(); it++)
            size += it->Multiplicity();
        return size;
    }

    void OutputSequences(string output_fname) const {
        ofstream out(output_fname.c_str());
        size_t id = 1;
        for(auto it = begin(); it != end(); it++) {
            out << ">antibody_" << id << endl;
            out << it->Sequence() << endl;
            id++;
        }
        out.close();
    }

    void OutputRepertoire(string output_fname) const {
        ofstream out(output_fname.c_str());
        size_t id = 1;
        for(auto it = begin(); it != end(); it++) {
            for(size_t i = 0; i < it->Multiplicity(); i++) {
                out << ">antibody_" << id << "_multiplicity_" << it->Multiplicity() << "_copy_" << i + 1 << endl;
                out << it->Sequence() << endl;
            }
            id++;
        }
        out.close();
    }

    void OutputVDJRecombination(string output_fname) const {
        ofstream out(output_fname.c_str());
        size_t id = 1;
        for(auto it = begin(); it != end(); it++) {
            out << "antibody_" << id << "\t" <<
                    it->IgVariableRegion()->VDJ_Recombination()->VJDRecombinationString() << endl;
            id++;
        }
        out.close();
    }

    void OutputMultiplicities(string output_fname) const {
        ofstream out(output_fname.c_str());
        for(auto it = begin(); it != end(); it++)
            out << it->Multiplicity() << endl;
        out.close();
    }

    void OutputSHMPositions(string output_fname) const {
        ofstream out(output_fname.c_str());
        for(auto it = begin(); it != end(); it++) {
            auto shm_settings = it->IgVariableRegion()->GetSHMSettings();
            size_t seq_length = it->Sequence().size();
            for(auto shm = shm_settings.begin(); shm != shm_settings.end(); shm++) {
                if(shm->first > seq_length)
                    continue;
                out << shm->first << "\t" << seq_length << endl;
            }

        }
    }
};

typedef Repertoire<HC_Cluster, HC_ClusterIterator> HC_Repertoire;
typedef Repertoire<LC_Cluster, LC_ClusterIterator> LC_Repertoire;

typedef shared_ptr<HC_Repertoire> HC_Repertoire_Ptr;
typedef shared_ptr<LC_Repertoire> LC_Repertoire_Ptr;