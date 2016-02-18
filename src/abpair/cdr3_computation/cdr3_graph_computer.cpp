#include <logger/logger.hpp>
#include <verify.hpp>
#include <queue>
#include "cdr3_graph_computer.hpp"
#include "simple_cdr3_calculator.hpp"

void Cdr3GraphComputer::Initialize() {
    for(size_t i = 0; i < pairing_data_storage_.size(); i++)
        if(pairing_data_storage_[i]->Complete())
            pairing_data_cdr3s_.push_back(PairingDataCdr3());
}

void Cdr3GraphComputer::ComputeHCdr3s() {
    size_t ind = 0;
    for (auto it = pairing_data_storage_.cbegin(); it != pairing_data_storage_.cend(); it++) {
        if (!(*it)->Complete())
            continue;
        auto hc_isotypes = (*it)->HcIsotypes();
        for (auto hc = hc_isotypes.begin(); hc != hc_isotypes.end(); hc++) {
            auto hc_seqs = (*it)->GetSequencesByIsotype(*hc);
            for (auto hc_seq = hc_seqs->cbegin(); hc_seq != hc_seqs->cend(); hc_seq++) {
                auto cdr3 = SimpleCdr3Calculator().FindCdr3Positions(*hc_seq);
                if (hcdr3_map_.find(cdr3) == hcdr3_map_.end())
                    hcdr3_map_[cdr3] = std::vector<size_t>();
                hcdr3_map_[cdr3].push_back(ind);
                pairing_data_cdr3s_[ind].hcdr3s.insert(cdr3);
            }
        }
        ind++;
    }
    INFO("HCDR3 map of size " << hcdr3_map_.size() << " was computed");
}

void Cdr3GraphComputer::ComputeKCdr3s() {
    size_t ind = 0;
    for(auto it = pairing_data_storage_.cbegin(); it != pairing_data_storage_.cend(); it++) {
        if (!(*it)->Complete())
            continue;
        if ((*it)->KappaChainCount() != 0) {
            auto kappa_seq = (*it)->GetSequencesByIsotype(IgIsotypeHelper::GetKappaIsotype());
            for (auto kit = kappa_seq->cbegin(); kit != kappa_seq->cend(); kit++) {
                auto cdr3 = SimpleCdr3Calculator().FindCdr3Positions(*kit);
                if (kcdr3_map_.find(cdr3) != kcdr3_map_.end())
                    kcdr3_map_[cdr3] = std::vector<size_t>();
                kcdr3_map_[cdr3].push_back(ind);
                pairing_data_cdr3s_[ind].kcdr3s.insert(cdr3);
            }
        }
        ind++;
    }
    INFO("KCDR3 map of size " << kcdr3_map_.size() << " was computed");
}

void Cdr3GraphComputer::ComputeLCdr3s() {
    size_t ind = 0;
    for(auto it = pairing_data_storage_.cbegin(); it != pairing_data_storage_.cend(); it++) {
        if (!(*it)->Complete())
            continue;
        if ((*it)->LambdaChainCount() != 0) {
            auto lambda_seq = (*it)->GetSequencesByIsotype(IgIsotypeHelper::GetLambdaIsotype());
            for (auto lit = lambda_seq->cbegin(); lit != lambda_seq->cend(); lit++) {
                auto cdr3 = SimpleCdr3Calculator().FindCdr3Positions(*lit);
                if (lcdr3_map_.find(cdr3) != lcdr3_map_.end())
                    lcdr3_map_[cdr3] = std::vector<size_t>();
                lcdr3_map_[cdr3].push_back(ind);
                pairing_data_cdr3s_[ind].lcdr3s.insert(cdr3);
            }
        }
        ind++;
    }
    INFO("LCDR3 map of size " << lcdr3_map_.size() << " was computed");
}

void Cdr3GraphComputer::InitializeMaps() {
    for(auto it = hcdr3_map_.begin(); it != hcdr3_map_.end(); it++)
        hcdr3_processed_[it->first] = false;
    for(auto it = kcdr3_map_.begin(); it != kcdr3_map_.end(); it++)
        kcdr3_processed_[it->first] = false;
    for(auto it = lcdr3_map_.begin(); it != lcdr3_map_.end(); it++)
        lcdr3_processed_[it->first] = false;
    for(size_t i = 0; i < pairing_data_cdr3s_.size(); i++)
        record_processed_.push_back(false);
}

std::string Cdr3GraphComputer::GetStartHCdr3() {
    for(auto it = hcdr3_processed_.begin(); it != hcdr3_processed_.end(); it++)
        if(!it->second)
            return it->first;
    return "";
}

std::string Cdr3GraphComputer::GetStartKCdr3() {
    for(auto it = kcdr3_processed_.begin(); it != kcdr3_processed_.end(); it++)
        if(!it->second)
            return it->first;
    return "";
}

std::string Cdr3GraphComputer::GetStartLCdr3() {
    for(auto it = lcdr3_processed_.begin(); it != lcdr3_processed_.end(); it++)
        if(!it->second)
            return it->first;
    return "";
}

void Cdr3GraphComputer::ComputeCdr3GraphFromStartCdr3(std::string start_hcdr3) {
    std::queue<std::string> hcdr3_queue;
    std::queue<std::string> kcdr3_queue;
    std::queue<std::string> lcdr3_queue;
    hcdr3_queue.push(start_hcdr3);
    std::string hcdr3;
    std::string kcdr3;
    std::string lcdr3;
    std::unordered_set<std::string> hcdr3s;
    std::unordered_set<std::string> kcdr3s;
    std::unordered_set<std::string> lcdr3s;
    INFO("Start HCDR3: " << start_hcdr3);
    while(!hcdr3_queue.empty() or !lcdr3_queue.empty() or !kcdr3_queue.empty()) {
        if (!hcdr3_queue.empty()) {
            hcdr3 = hcdr3_queue.front();
            hcdr3_queue.pop();
        }
        if (!kcdr3_queue.empty()) {
            kcdr3 = kcdr3_queue.front();
            kcdr3_queue.pop();
        }
        if (!lcdr3_queue.empty()) {
            lcdr3 = lcdr3_queue.front();
            lcdr3_queue.pop();
        }
        std::vector<size_t> indices;
        if (hcdr3 != "") {
            INFO("Current H CDR3: " << hcdr3);
            hcdr3_processed_[hcdr3] = true;
            hcdr3s.insert(hcdr3);
            indices.insert(indices.end(), hcdr3_map_[hcdr3].begin(), hcdr3_map_[hcdr3].end());
        }
        if (kcdr3 != "") {
            kcdr3_processed_[kcdr3] = true;
            kcdr3s.insert(kcdr3);
            indices = kcdr3_map_[kcdr3];
            indices.insert(indices.end(), kcdr3_map_[hcdr3].begin(), kcdr3_map_[hcdr3].end());
        }
        if (lcdr3 != "") {
            lcdr3_processed_[lcdr3] = true;
            lcdr3s.insert(lcdr3);
            indices = lcdr3_map_[lcdr3];
            indices.insert(indices.end(), lcdr3_map_[hcdr3].begin(), lcdr3_map_[hcdr3].end());
        }
        INFO("# indices: " << indices.size());
        for (auto ind = indices.begin(); ind != indices.end(); ind++) {
            // add other hc cdr3s
            auto hcdr3_vector = pairing_data_cdr3s_[*ind].hcdr3s;
            for (auto it = hcdr3_vector.begin(); it != hcdr3_vector.end(); it++)
                if (!hcdr3_processed_[*it]) {
                    hcdr3_processed_[*it] = true;
                    hcdr3_queue.push(*it);
                    INFO("H CDR3 " << *it << "-> set");
                }
            // add other kappa cdr3s
            auto kcdr3_vector = pairing_data_cdr3s_[*ind].kcdr3s;
            for (auto it = kcdr3_vector.begin(); it != kcdr3_vector.end(); it++)
                if (!kcdr3_processed_[*it]) {
                    kcdr3_processed_[*it] = true;
                    kcdr3_queue.push(*it);
                    INFO("K CDR3 " << *it << "-> set");
                }
            // add other lambda cdr3s
            auto lcdr3_vector = pairing_data_cdr3s_[*ind].lcdr3s;
            for (auto it = lcdr3_vector.begin(); it != lcdr3_vector.end(); it++)
                if (!lcdr3_processed_[*it]) {
                    lcdr3_processed_[*it] = true;
                    lcdr3_queue.push(*it);
                    INFO("L CDR3 " << *it << "-> set");
                }
        }
    }
    INFO("# HCDR3: " << hcdr3s.size() << ", # KCDR3: " << kcdr3s.size() << ", # LCDR3: " << lcdr3.size());
    related_hcdr3_.push_back(hcdr3s);
    related_kcdr3_.push_back(kcdr3s);
    related_lcdr3_.push_back(lcdr3s);
}

void Cdr3GraphComputer::ComputeCdr3Graphs() {
    InitializeMaps();
    std::string hcdr3 = GetStartHCdr3();
    while(hcdr3 != "") {
        ComputeCdr3GraphFromStartCdr3(hcdr3);
        break;
        hcdr3 = GetStartHCdr3();
    }
}

void Cdr3GraphComputer::Compute() {
    ComputeHCdr3s();
    ComputeKCdr3s();
    ComputeLCdr3s();
    ComputeCdr3Graphs();
}
