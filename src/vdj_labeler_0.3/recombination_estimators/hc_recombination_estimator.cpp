#include "logger/logger.hpp"
#include "hc_recombination_estimator.hpp"

void HcRecombinationEstimator::Update(HcRecombinationStoragePtr recombination_storage) {
    num_recombinations_.push_back(recombination_storage->size());
    if(recombination_storage->size() == 0)
        return;
    size_t min_num_shms = size_t(-1);
    size_t min_shms_index = size_t(-1);
    for(size_t i = 0; i < recombination_storage->size(); i++)
        if((*recombination_storage)[i].SHMsNumber() < min_num_shms) {
            min_num_shms = (*recombination_storage)[i].SHMsNumber();
            min_shms_index = i;
        }
    min_number_shms_.push_back(min_num_shms);
    v_event_lens_.push_back((*recombination_storage)[min_shms_index].V().RightCleavageLength());
    d_levent_lens_.push_back((*recombination_storage)[min_shms_index].D().LeftCleavageLength());
    d_revent_lens_.push_back((*recombination_storage)[min_shms_index].D().RightCleavageLength());
    j_event_lens_.push_back((*recombination_storage)[min_shms_index].J().LeftCleavageLength());
}

void HcRecombinationEstimator::OutputRecombinationNumber() {
    // todo: refactor it and put files in dir
    std::string fname = "recombination_number_dist.txt";
    std::ofstream out(fname);
    for(auto it = num_recombinations_.begin(); it != num_recombinations_.end(); it++)
        out << *it << std::endl;
    out.close();
    INFO("Distribution of recombination numbers was written to " << fname);
}

void HcRecombinationEstimator::OutputSHMsDistribution() {
    // todo: refactor it and put files in dir
    std::string fname = "shms_number_dist.txt";
    std::ofstream out(fname);
    for(auto it = min_number_shms_.begin(); it != min_number_shms_.end(); it++)
        out << *it << std::endl;
    out.close();
    INFO("Distribution of # min SHMs was written to " << fname);
}

void HcRecombinationEstimator::OutputRecombinationEvents() {
    std::string fname = "min_shms_event_dist.txt";
    std::ofstream out(fname);
    for(size_t i = 0; i < v_event_lens_.size(); i++)
        out << v_event_lens_[i] << "\t" << d_levent_lens_[i] << "\t" << d_revent_lens_[i] <<
                "\t" << j_event_lens_[i] << std::endl;
    out.close();
    INFO("Distribution of recombination event lengths for min SHMs was written to " << fname);
}