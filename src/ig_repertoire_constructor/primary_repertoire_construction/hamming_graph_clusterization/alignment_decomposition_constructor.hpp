#pragma once

#include "hg_decomposition.hpp"

class PositionCounter {
    size_t a_;
    size_t c_;
    size_t g_;
    size_t t_;

    vector<size_t> a_indices_;
    vector<size_t> c_indices_;
    vector<size_t> g_indices_;
    vector<size_t> t_indices_;

    void UpdateA(size_t index, size_t multiplicity = 1) {
        a_ += multiplicity;
        a_indices_.push_back(index);
    }

    void UpdateC(size_t index, size_t multiplicity = 1) {
        c_ += multiplicity;
        c_indices_.push_back(index);
    }

    void UpdateG(size_t index, size_t multiplicity = 1) {
        g_ += multiplicity;
        g_indices_.push_back(index);
    }

    void UpdateT(size_t index, size_t multiplicity = 1) {
        t_ += multiplicity;
        t_indices_.push_back(index);
    }

public:
    PositionCounter() :
        a_(0), c_(0), g_(0), t_(0) { }

    void Update(char nucl, size_t index, size_t multiplicity = 1) {
        if(nucl == 'A')
            UpdateA(index, multiplicity);
        else if(nucl == 'C')
            UpdateC(index, multiplicity);
        else if(nucl == 'G')
            UpdateG(index, multiplicity);
        else if(nucl == 'T')
            UpdateT(index, multiplicity);
    }

    size_t A() const { return a_; }
    size_t C() const { return c_; }
    size_t G() const { return g_; }
    size_t T() const { return t_; }

    size_t Sum() const { return A() + C() + G() + T(); }

    double PercA() const { return double(A()) / double(Sum()); }
    double PercC() const { return double(C()) / double(Sum()); }
    double PercG() const { return double(G()) / double(Sum()); }
    double PercT() const { return double(T()) / double(Sum()); }

    bool IsMismatch() const { return A() * C() != 0 or A() * G() != 0 or A() * T() != 0 or C() * G() != 0 or
                C() * T() != 0 or G() * T() != 0; }

    bool MismatchIsTrivial() const {
        return (A() == Sum() - 1 or C() == Sum() - 1 or G() == Sum() - 1 or T() == Sum() - 1);
    }

    const vector<size_t>& IndicesA() const { return a_indices_; }
    const vector<size_t>& IndicesC() const { return c_indices_; }
    const vector<size_t>& IndicesG() const { return g_indices_; }
    const vector<size_t>& IndicesT() const { return t_indices_; }

    size_t NonzeroPerc() const {
        size_t res = 0;
        if(A() != 0) res++;
        if(C() != 0) res++;
        if(G() != 0) res++;
        if(T() != 0) res++;
        return res;
    }

    double MaxPerc() const {
        return max<double>(max<double>(PercA(), PercC()), max<double>(PercG(), PercT()));
    }

    char MaxNucl() const {
        double epsilon = 0.0001;
        if(fabs(MaxPerc() - PercA()) < epsilon)
            return 'A';
        if(fabs(MaxPerc() - PercC()) < epsilon)
            return 'C';
        if(fabs(MaxPerc() - PercG()) < epsilon)
            return 'G';
        return 'T';
    }

    double SecondMaxPerc() const {
        char max_nucl = MaxNucl();
        double second_max = 0;
        if(PercA() >= second_max and max_nucl != 'A')
            second_max = PercA();
        if(PercC() >= second_max and max_nucl != 'C')
            second_max = PercC();
        if(PercG() >= second_max and max_nucl != 'G')
            second_max = PercG();
        if(PercT() >= second_max and max_nucl != 'T')
            second_max = PercT();
        return second_max;
    }

};

ostream& operator<<(ostream &out, const PositionCounter &counter) {
    out << "A: " << counter.A() << " (" << counter.PercA() << "%), C: " << counter.C() << " (" <<
            counter.PercC() << "%), G: " << counter.G() << " (" << counter.PercG() << "%), T: " <<
            counter.T() << " (" << counter.PercT() << "%)";
    out << endl << "A indices: ";
    for(auto it = counter.IndicesA().begin(); it != counter.IndicesA().end(); it++)
        out << *it << " ";
    out << endl << "C indices: ";
    for(auto it = counter.IndicesC().begin(); it != counter.IndicesC().end(); it++)
        out << *it << " ";
    out << endl << "G indices: ";
    for(auto it = counter.IndicesG().begin(); it != counter.IndicesG().end(); it++)
        out << *it << " ";
    out << endl << "T indices: ";
    for(auto it = counter.IndicesT().begin(); it != counter.IndicesT().end(); it++)
        out << *it << " ";
    return out;
}

class PositionCountersDecomposer {
    // input parameters
    const vector<PositionCounter> &counters_;
    HG_DecompositionPtr decomposition_;
    double min_mismatch_perc_;

    // auxiliary structs
    vector<size_t> valid_counter_indices_;
    vector<double> counter_scores_;


    bool AllPercLarge(const PositionCounter &counter) {
        if(counter.A() != 0 and counter.PercA() < min_mismatch_perc_)
            return false;
        if(counter.C() != 0 and counter.PercC() < min_mismatch_perc_)
            return false;
        if(counter.G() != 0 and counter.PercG() < min_mismatch_perc_)
            return false;
        if(counter.T() != 0 and counter.PercT() < min_mismatch_perc_)
            return false;
        return true;
    }

    bool AllPercLargeAndOneSmall(const PositionCounter &counter) {
        if(counter.NonzeroPerc() < 3)
            return false;
        size_t num_large_perc = 0;
        if(counter.PercA() >= min_mismatch_perc_)
            num_large_perc++;
        if(counter.PercC() >= min_mismatch_perc_)
            num_large_perc++;
        if(counter.PercG() >= min_mismatch_perc_)
            num_large_perc++;
        if(counter.PercT() >= min_mismatch_perc_)
            num_large_perc++;
        return num_large_perc >= 2;
    }

    bool CounterIsValid(const PositionCounter &counter) {
        if(AllPercLarge(counter))
            return true;
        return AllPercLargeAndOneSmall(counter);

    }

    double ComputeCounterScore(const PositionCounter &counter) {
        //TRACE("Maximal %: " << counter.MaxPerc() << "\t" << counter.SecondMaxPerc());
        return 1 - (counter.MaxPerc() - counter.SecondMaxPerc());
    }

    void ComputeCounterScores() {
        for(size_t i = 0; i < counters_.size(); i++)
            if(CounterIsValid(counters_[i])) {
                valid_counter_indices_.push_back(i);
                double score = ComputeCounterScore(counters_[i]);
                counter_scores_.push_back(score);
                TRACE("Position counter #" << i << " is valid, score: " << score);
            }
    }

    size_t GetBestCounter() {
        double score = 0;
        size_t index = size_t(-1);
        for(size_t i = 0; i < counter_scores_.size(); i++)
            if(counter_scores_[i] > score) {
                score = counter_scores_[i];
                index = valid_counter_indices_[i];
            }
        return index;
    }

    void CreateTrivialDecomposition() {
        size_t cur_class_id = decomposition_->Size();
        for(auto it = counters_[0].IndicesA().begin(); it != counters_[0].IndicesA().end(); it++)
            decomposition_->SetClass(*it, cur_class_id);
        for(auto it = counters_[0].IndicesC().begin(); it != counters_[0].IndicesC().end(); it++)
            decomposition_->SetClass(*it, cur_class_id);
        for(auto it = counters_[0].IndicesG().begin(); it != counters_[0].IndicesG().end(); it++)
            decomposition_->SetClass(*it, cur_class_id);
        for(auto it = counters_[0].IndicesT().begin(); it != counters_[0].IndicesT().end(); it++)
            decomposition_->SetClass(*it, cur_class_id);
    }

    void CreateDecompositionByIndex(size_t index) {
        size_t cur_class_id = decomposition_->Size();
        if(counters_[index].A() != 0) {
            for(auto it = counters_[index].IndicesA().begin(); it != counters_[index].IndicesA().end(); it++)
                decomposition_->SetClass(*it, cur_class_id);
            cur_class_id++;
        }
        if(counters_[index].C() != 0) {
            for(auto it = counters_[index].IndicesC().begin(); it != counters_[index].IndicesC().end(); it++)
                decomposition_->SetClass(*it, cur_class_id);
            cur_class_id++;
        }
        if(counters_[index].G() != 0) {
            for(auto it = counters_[index].IndicesG().begin(); it != counters_[index].IndicesG().end(); it++)
                decomposition_->SetClass(*it, cur_class_id);
            cur_class_id++;
        }
        if(counters_[index].T() != 0) {
            for(auto it = counters_[index].IndicesT().begin(); it != counters_[index].IndicesT().end(); it++)
                decomposition_->SetClass(*it, cur_class_id);
            cur_class_id++;
        }
    }

public:
    PositionCountersDecomposer(const vector<PositionCounter> &counters,
            HG_DecompositionPtr decomposition,
            double min_mismatch_perc) :
            counters_(counters),
            decomposition_(decomposition),
            min_mismatch_perc_(min_mismatch_perc) {}

    HG_DecompositionPtr CreateDecomposition() {
        TRACE(counters_.size() << " non trivial counters were detected");
        ComputeCounterScores();
        TRACE("# significant counters: " << valid_counter_indices_.size());
        if(valid_counter_indices_.size() == 0) {
            TRACE("Trivial secomposition will be returned");
            CreateTrivialDecomposition();
        }
        else {
            size_t best_mismatch = GetBestCounter();
            assert(best_mismatch != size_t(-1));
            TRACE("Index of best mismatch is " << best_mismatch);
            CreateDecompositionByIndex(best_mismatch);
        }
        return decomposition_;
    }

private:
    DECL_LOGGER("PositionCountersDecomposer");
};

class AlignmentDecompositionConstructor {
    // input parameters
    CRS_HammingGraph_Ptr hamming_graph_ptr_;
    HG_CollapsedStructs_Ptr collapsed_struct_ptr_;
    HG_DecompositionPtr primary_decomposition_;
    const vector<ig_repertoire_constructor::SplicedRead> &reads_;
    size_t group_id_;

    // output decomposition
    HG_DecompositionPtr output_decomposition_;

    PositionCounter CountPosition(const set<size_t> &cur_class, const vector<size_t> &read_indices, size_t pos) {
        PositionCounter mismatch_counter;
        size_t read_index = 0;
        for(auto it = cur_class.begin(); it != cur_class.end(); it++) {
            mismatch_counter.Update(nucl(reads_[read_indices[read_index]][pos]), *it,
                    collapsed_struct_ptr_->MultiplicityOfNewVertex(*it));
            read_index++;
        }
        return mismatch_counter;
    }

    void CreateTrivialDecomposition(const set<size_t> &cur_class) {
        size_t cur_class_id = output_decomposition_->Size();
        for(auto it = cur_class.begin(); it != cur_class.end(); it++) {
            //TRACE("Cur vertex: " << *it << ", # vertices: " << output_decomposition_->VertexNumber());
            output_decomposition_->SetClass(*it, cur_class_id);
        }
    }

    void CreateClassDecomposition(const vector<PositionCounter> &counters) {
        double min_mismatch_perc = .2; // todo: compute and move to config
        PositionCountersDecomposer decomposer(counters, output_decomposition_, min_mismatch_perc);
        output_decomposition_ = decomposer.CreateDecomposition();
    }

    void ProcessClass(const set<size_t> &cur_class, size_t class_id) {
        vector<size_t> old_read_ids;
        for(auto it = cur_class.begin(); it != cur_class.end(); it++)
            old_read_ids.push_back(collapsed_struct_ptr_->OldVerticesList()[*it]);
        size_t read_size = reads_[0].size();
        size_t num_trivial_mismatches = 0;
        vector<PositionCounter> non_trivial_counters;
        // tmp
        ofstream out("alignment_decomposition.txt", std::ios_base::app);
        for(size_t i = 0; i < read_size; i++) {
            PositionCounter counter = CountPosition(cur_class, old_read_ids, i);
            if(counter.IsMismatch()) {
                if(counter.MismatchIsTrivial()) {
                    num_trivial_mismatches++;
                }
                else {
                    TRACE("Pos: " << i << ", non trivial mismatch: " << counter);
                    out << group_id_ << "\t" << class_id << "\t" <<
                            counter.A() << "\t" << counter.C() << "\t" <<
                            counter.G() << "\t" << counter.T() << "\t" <<
                            primary_decomposition_->RealSizeOfClass(class_id, collapsed_struct_ptr_) <<
                            "\t" << i << "\t" << read_size << endl;
                    non_trivial_counters.push_back(counter);
                }
            }
        }
        if(non_trivial_counters.size() == 0)
            CreateTrivialDecomposition(cur_class);
        else
            CreateClassDecomposition(non_trivial_counters);
        TRACE("# trivial mismatches: " << num_trivial_mismatches);
    }

public:
    AlignmentDecompositionConstructor(CRS_HammingGraph_Ptr hamming_graph_ptr,
            HG_CollapsedStructs_Ptr collapsed_struct_ptr,
            HG_DecompositionPtr primary_decomposition,
            const vector<ig_repertoire_constructor::SplicedRead> &reads,
            size_t group_id) :
        hamming_graph_ptr_(hamming_graph_ptr),
        collapsed_struct_ptr_(collapsed_struct_ptr),
        primary_decomposition_(primary_decomposition),
        reads_(reads),
        group_id_(group_id),
        output_decomposition_(new HG_Decomposition(primary_decomposition->VertexNumber())) { }

    HG_DecompositionPtr ConstructDecomposition() {
        size_t min_class_size = 2; // todo: define and move to config
        for(size_t i = 0; i < primary_decomposition_->Size(); i++) {
            auto cur_class = primary_decomposition_->GetClass(i);
            TRACE("------------");
            TRACE("Processing class #" << i << ", size: " << cur_class.size());
            if(cur_class.size() < min_class_size) {
                TRACE("Class is trivial, decomposition will be skipped");
                CreateTrivialDecomposition(cur_class);
                continue;
            }
            ProcessClass(cur_class, i);
        }
        TRACE("# constructed classes: " << output_decomposition_->Size());
        return output_decomposition_;
    }

private:
    DECL_LOGGER("AlignmentDecompositionConstructor");
};