#pragma once

#include <set>
#include <omp.h>

#include "shared_kmer_string_set.hpp"
#include "repertoire.hpp"
#include "logger/logger.hpp"

namespace ig_repertoire_constructor {

using string_inexact_match_index::AllIndexSelector;
using string_inexact_match_index::SharedKMerString;
using string_inexact_match_index::SharedKMerStringSet;

class SingletonGluer {
    DECL_LOGGER("SingletonGluer");

    struct IndexedShift {
        size_t cluster_ind;
        int relative_shift;
    };

    static size_t max_kmer_occurences_;
    static size_t k_;
    static size_t max_dist_;
    static size_t min_overlap_;
    static size_t common_mismatches_threshold_;
    static size_t diverse_mutations_threshold_;

    Repertoire & repertoire_;
    set<size_t> const & bad_clusters_ids_;

    map<size_t, IndexedShift> gluing_map_;
    map<size_t, size_t> distance_hist_;

    set<size_t> successfully_glued_;
    size_t candidate_to_gluing_;

    omp_lock_t lock_;

    bool ClusterShouldBeGlued(size_t index) {
        return (repertoire_.GetClusterByIndex(index).size() == 1 ||
                bad_clusters_ids_.find(repertoire_.GetClusterByIndex(index).id()) != bad_clusters_ids_.end());
    }

    bool GluingCandidateGood(size_t index) {
        return (!ClusterShouldBeGlued(index));
    }

    void PrintGluingMapTraces() const;

    void CreateGluingMap();
    void UpdateGluingMap(SharedKMerStringSet const &index, vector <size_t> const &index_seqnum_to_cluster_ind,
                         size_t to_glue);

    bool CheckMismatches(size_t max_mismatches, std::vector <std::vector <std::pair <size_t, char> > > const &mismatches_for_singleton) const;
    bool TryToGlue(size_t gluing_dst, const vector<size_t> & gluing_srcs);

    void UpdateRepertoire();

    /*
    void OutputStats() {
        ofstream out(fnames_config_.output_stats.c_str()); 
        for(auto it = distance_hist_.begin(); it != distance_hist_.end(); it++)
            for(size_t i = 0; i < it->second; i++)
                out << it->first << endl;
        cout << "Statistics for glued singletons were written to " << fnames_config_.output_stats << endl;
        out.close();
    }
    */

    size_t GetShiftHammingDistance(size_t relative_shift, const std::string & s1, const std::string & s2);

public:
    SingletonGluer(RepertoirePtr repertoire_ptr,
        std::set <size_t> const & bad_clusters_ids) :
        repertoire_(*repertoire_ptr),
        bad_clusters_ids_(bad_clusters_ids),
        candidate_to_gluing_(0) { }

    static void init_singleton_gluer(size_t k, size_t max_kmer_occurences,
                                     size_t max_distance, size_t min_overlap,
                                     size_t common_mismatches_threshold, size_t diverse_mutations_threshold) {
        k_ = k;
        max_kmer_occurences_ = max_kmer_occurences;
        max_dist_ = max_distance;
        min_overlap_ = min_overlap;
        common_mismatches_threshold_ = common_mismatches_threshold;
        diverse_mutations_threshold_ = diverse_mutations_threshold;
    }

    void GlueSingletons() { 
        CreateGluingMap();
        UpdateRepertoire();
        // OutputStats();
    }
};

size_t SingletonGluer::max_kmer_occurences_ = 5000;
size_t SingletonGluer::k_ = 60;
size_t SingletonGluer::max_dist_ = 4;
size_t SingletonGluer::min_overlap_ = 200;
size_t SingletonGluer::common_mismatches_threshold_ = 2;
size_t SingletonGluer::diverse_mutations_threshold_ = 1;

void SingletonGluer::PrintGluingMapTraces() const {
    TRACE("Gluing map:");
    TRACE("ID\t<-\tID");
    for(auto it = gluing_map_.begin(); it != gluing_map_.end(); it++) {
        TRACE(repertoire_.GetClusterByIndex(it->first).id() << "\t\t" <<
              repertoire_.GetClusterByIndex(it->second.cluster_ind).id());
    }

    TRACE("Distance hist:");
    TRACE("Dist\tNum");
    for(auto it = distance_hist_.begin(); it != distance_hist_.end(); it++)
        TRACE(it->first << "\t" << it->second);
}

void SingletonGluer::CreateGluingMap() {
    std::shared_ptr <AllIndexSelector> selector(new AllIndexSelector());
    SharedKMerStringSet index(selector, k_, min_overlap_, max_kmer_occurences_);
    vector <size_t> index_seqnum_to_cluster_ind;
    for (size_t i = 0; i != repertoire_.size(); ++i) {
        if (GluingCandidateGood(i)) {
            index.Insert(repertoire_.GetClusterByIndex(i).sequence().str());
            index_seqnum_to_cluster_ind.push_back(i);
        }
    }

    omp_init_lock(&lock_);
    #pragma omp parallel for num_threads(ig_cfg::get().rp.threads_count)
    for(size_t i = 0; i < repertoire_.size(); i++) {
        if(ClusterShouldBeGlued(i)) {
            candidate_to_gluing_++;
            UpdateGluingMap(index, index_seqnum_to_cluster_ind, i);
        }
    }
    omp_destroy_lock(&lock_);

    PrintGluingMapTraces();

    INFO(successfully_glued_.size() << " clusters from " << candidate_to_gluing_ << " can be glued");
}

void SingletonGluer::UpdateGluingMap(SharedKMerStringSet const &index, vector <size_t> const &index_seqnum_to_cluster_ind,
                                     size_t to_glue) {
    vector <SharedKMerString> candidates;
    index.Find(repertoire_.GetClusterByIndex(to_glue).sequence().str(), candidates);

    size_t best_dist = max_dist_ + 1;
    IndexedShift best_candidate = {0, 0};
    bool candidate_found = false;
    for(auto it = candidates.begin(); it != candidates.end(); it++) {
        size_t candidate_ind = index_seqnum_to_cluster_ind[it->seq_number_in_set];
        if(GluingCandidateGood(candidate_ind) and candidate_ind != to_glue) {
            size_t dist = GetShiftHammingDistance(it->pattern_start_position,
                                    repertoire_.GetClusterByIndex(candidate_ind).sequence().str(),
                                    repertoire_.GetClusterByIndex(to_glue).sequence().str());

            if(dist < best_dist) {
                candidate_found = true;
                best_dist = dist;
                best_candidate = {candidate_ind, it->pattern_start_position};
            }
        }
    }
    if(candidate_found) {
        omp_set_lock(&lock_);
        gluing_map_[to_glue] = best_candidate;
        distance_hist_[best_dist]++;
        successfully_glued_.insert(repertoire_.GetClusterByIndex(to_glue).id());
        omp_unset_lock(&lock_);
        TRACE("Shift " << best_candidate.relative_shift);
    }
}

bool SingletonGluer::CheckMismatches(size_t max_mismatches, const std::vector<std::vector<std::pair<size_t, char> > > &mismatches_for_singleton) const {
    size_t common_mismatches = 0;
    for (auto mismch_it = mismatches_for_singleton[0].begin();
            mismch_it != mismatches_for_singleton[0].end(); ++mismch_it) {
        for (size_t singleton_i = 1; singleton_i != mismatches_for_singleton.size();
                ++singleton_i) {
            if (std::find(mismatches_for_singleton[singleton_i].begin(),
                          mismatches_for_singleton[singleton_i].end(), *mismch_it) !=
                    mismatches_for_singleton[singleton_i].end()) {
                ++common_mismatches;
            }
            /*
            std::ostringstream str_stream;
            str_stream << *lpos_it << " ";
            TRACE(str_stream.str());
            */
        }
    }

    if (max_mismatches - common_mismatches <= diverse_mutations_threshold_ &&
            common_mismatches >= common_mismatches_threshold_) {
        return false;
    }
    return true;
}

bool SingletonGluer::TryToGlue(size_t gluing_dst, const vector<size_t> & gluing_srcs) {
    //TODO: need to refactor
    std::string cluster_seq = repertoire_.GetClusterByIndex(gluing_dst).sequence().str();
    int start_nucl_added = 0;
    TRACE("Gluing " << repertoire_.GetClusterByIndex(gluing_dst).id() << " cluster");
    std::vector <std::vector <std::pair <size_t, char> > > mismatches_for_singleton;
    size_t max_mismatches = 0;
    for (auto src_it = gluing_srcs.begin(); src_it != gluing_srcs.end(); ++src_it) {
        int curr_shift = gluing_map_[*src_it].relative_shift + start_nucl_added;
        std::string glued_seq = repertoire_.GetClusterByIndex(*src_it).sequence().str();
        if (GetShiftHammingDistance(curr_shift, cluster_seq, glued_seq) > max_dist_) {
            return false;
        }

        std::vector <std::pair <size_t, char> > curr_mismatches;

        size_t pos1 = (curr_shift > 0) ? curr_shift: 0;
        size_t pos2 = (curr_shift > 0) ? 0 : -curr_shift;
        size_t overlap_size = min(cluster_seq.size() - pos1, glued_seq.size() - pos2);
        for (size_t i = 0; i != overlap_size; ++i) {
            if (cluster_seq[pos1 + i] != glued_seq[pos2 + i]) {
                curr_mismatches.push_back(std::make_pair(pos1 + i, glued_seq[pos2 + i]));
            }
        }

        if (curr_mismatches.size() > max_mismatches) {
            max_mismatches = curr_mismatches.size();
        }
        mismatches_for_singleton.push_back(curr_mismatches);

        if (curr_shift < 0) {
            cluster_seq.insert(0, glued_seq, 0, abs(curr_shift));
            start_nucl_added += abs(curr_shift);
            for (auto singleton_it = mismatches_for_singleton.begin();
                    singleton_it != mismatches_for_singleton.end(); ++singleton_it) {
                for (auto it = singleton_it->begin(); it != singleton_it->end(); ++it) {
                    it->first -= abs(curr_shift);
                }
            }
        }
        int glue_after_seq = int(glued_seq.size() - cluster_seq.size());
        if (glue_after_seq > 0) {
            cluster_seq += glued_seq.substr(glued_seq.size() - glue_after_seq);
        }
    }

    if (!CheckMismatches(max_mismatches, mismatches_for_singleton)) {
        return false;
    }

    if (gluing_srcs.size() != 1) {
        TRACE("Glued " << gluing_srcs.size());
    }

    int new_reading_frame = repertoire_.GetClusterByIndex(gluing_dst).reading_frame();
    if (new_reading_frame >= 0) {
        new_reading_frame = (new_reading_frame + start_nucl_added) % 3;
    }
    repertoire_.MergeClusters(gluing_dst, gluing_srcs, cluster_seq, new_reading_frame);
    return true;
}

void SingletonGluer::UpdateRepertoire() {
    map<size_t, vector<size_t> > clusters_to_glue;
    for (auto it = gluing_map_.begin(); it != gluing_map_.end(); it++) {
        size_t gluing_src = it->first;
        size_t gluing_dst = it->second.cluster_ind;
        clusters_to_glue[gluing_dst].push_back(gluing_src);
    }

    for (auto dst_it = clusters_to_glue.begin(); dst_it != clusters_to_glue.end(); ++dst_it) {
        if (!TryToGlue(dst_it->first, dst_it->second)) {
            for (auto src_it = dst_it->second.begin(); src_it != dst_it->second.end(); ++src_it) {
                successfully_glued_.erase(repertoire_.GetClusterByIndex(*src_it).id());
            }
        }
    }
    repertoire_.FilterIds(successfully_glued_);
    INFO(successfully_glued_.size() << " clusters from " << candidate_to_gluing_ << " were successfully glued");
}

size_t SingletonGluer::GetShiftHammingDistance(size_t relative_shift, const std::string & s1, const std::string & s2) {
    size_t pos1 = (relative_shift > 0) ? relative_shift : 0;
    size_t pos2 = (relative_shift > 0) ? 0 : -relative_shift;
    size_t overlap_size = min(s1.size() - pos1, s2.size() - pos2);
    size_t distance = 0;
    if (overlap_size < min_overlap_) {
        return max_dist_ + 1;
    }
    for (size_t i = 0; i != overlap_size; ++i) {
        if (s1[pos1 + i] != s2[pos2 + i]) {
            ++distance;
        }
        if (distance > max_dist_) {
            return distance;
        }
    }
    return distance;
}

}
