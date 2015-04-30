#pragma once

#include <repertoire.hpp>

namespace ig_repertoire_constructor {

struct stop_codon_pos {
    bool in_all_orfs;
    int predicted_reading_frame;
    size_t pos1;
    size_t pos2;
    size_t pos3;

    void SetIsInAllORFs(size_t seq_size) { in_all_orfs = pos1 < seq_size && pos2 < seq_size && pos3 < seq_size; }

    void FindReadingFrame(size_t seq_size) {
        if (pos1 == seq_size && pos2 < seq_size && pos3 < seq_size) {
            predicted_reading_frame = 0;
        } else if (pos2 == seq_size && pos1 < seq_size && pos3 < seq_size) {
            predicted_reading_frame = 1;
        } else if (pos3 == seq_size && pos1 < seq_size && pos2 < seq_size) {
            predicted_reading_frame = 2;
        } else {
            predicted_reading_frame = -1;
        }
    }

    bool operator==(const stop_codon_pos& obj) {
        return pos1 == obj.pos1 && pos2 == obj.pos2 && pos3 && obj.pos3;
    }
};

class StopCodonSearcher {
    static size_t StopCodonOccurs(string const & seq, size_t shift) {
        string stop_codon1 = "TAG";
        string stop_codon2 = "TAA";
        string stop_codon3 = "TGA";

        size_t start_pos = 3;
        size_t end_pos = seq.size() - 3;

        size_t start_aa = start_pos / 3 + (start_pos % 3 != 0) ? 1 : 0;
        size_t num_aa = (end_pos - shift) / 3;

        for(size_t i = start_aa; i < num_aa; i++) {
            size_t pos = i * 3 + shift;
            string aa = seq.substr(pos, 3);
            if(aa == stop_codon1 || aa == stop_codon2 || aa == stop_codon3)
                return pos;
        }
        return seq.size();
    }

public:
    static stop_codon_pos SequenceContainsStopCodon(string const & seq) {
        stop_codon_pos res;
        res.pos1 = StopCodonOccurs(seq, 0);
        res.pos2 = StopCodonOccurs(seq, 1);
        res.pos3 = StopCodonOccurs(seq, 2);
        res.SetIsInAllORFs(seq.size());
        res.FindReadingFrame(seq.size());
        return res;
    }
};

class AAVerificator {
    DECL_LOGGER("AAVerificator");
    std::set<size_t> incorrect_clusters_ids_;
    // vector<size_t> incorrect_clusters_sizes_;

public:
    void VerifyAndPredictReadingFrame(RepertoirePtr repertoire_ptr) {
        size_t not_identified_rf = 0;
        for(auto it = repertoire_ptr->begin(); it != repertoire_ptr->end(); it++) {
            stop_codon_pos pos = StopCodonSearcher::SequenceContainsStopCodon(it->sequence().str());
            if(pos.in_all_orfs) {
                TRACE("stop codon positions(" << it->sequence().size() << "): " <<
                      pos.pos1 << ", " << pos.pos2 << ", " << pos.pos3);
                incorrect_clusters_ids_.insert(it->id());
                // incorrect_clusters_sizes_.push_back(it->size);
            } else if (pos.predicted_reading_frame < 0) {
                not_identified_rf += 1;
            }
            it->reading_frame() = pos.predicted_reading_frame;
        }
        INFO(incorrect_clusters_ids_.size() << " clusters contain stop codons");
        INFO("Cannot identify reading frame for " << not_identified_rf << " clusters");
    }

    std::set <size_t> const & GetIncorrectClusterIds() const {
        return incorrect_clusters_ids_;
    }
};

}
