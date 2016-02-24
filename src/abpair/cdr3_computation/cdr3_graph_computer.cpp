#include <logger/logger.hpp>
#include <verify.hpp>
#include <queue>
#include "cdr3_graph_computer.hpp"
#include "simple_cdr3_calculator.hpp"

std::string GetHeaderForUmiSequence(const DropletBarcode &db, const IsotypeUmiSequence &umi_sequence) {
    std::stringstream ss;
    ss << "DB:" << db << "|UMI:" << umi_sequence.umi <<
    "|ISOTYPE:" << umi_sequence.isotype.str() << "|COUNT:" << umi_sequence.size;
    return ss.str();
}

size_t AnnotationHasher::operator()(const Annotation &obj) const {
    return std::hash<std::string>()(obj.cdr3) * std::hash<std::string>()(obj.v) * std::hash<std::string>()(obj.j);
}

void Cdr3GraphComputer::Initialize() {
    for(size_t i = 0; i < pairing_data_storage_.size(); i++)
        if(pairing_data_storage_[i]->Complete()) {
            pairing_data_cdr3s_.push_back(PairingDataCdr3());
            complete_indices_.push_back(i);
        }
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
                auto cdr3_seq = SimpleCdr3Calculator().FindCdr3Positions(*hc_seq);

                Annotation cdr3;
                cdr3.cdr3 = cdr3_seq; cdr3.v = hc_seq->v_gene; cdr3.j = hc_seq->j_gene;

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
                auto cdr3_seq = SimpleCdr3Calculator().FindCdr3Positions(*kit);

                Annotation cdr3;
                cdr3.cdr3 = cdr3_seq; cdr3.v = kit->v_gene; cdr3.j = kit->j_gene;

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
                auto cdr3_seq = SimpleCdr3Calculator().FindCdr3Positions(*lit);

                Annotation cdr3;
                cdr3.cdr3 = cdr3_seq; cdr3.v = lit->v_gene; cdr3.j = lit->j_gene;

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

size_t Cdr3GraphComputer::GetStartRecordIndex() {
    for(size_t i = 0; i != record_processed_.size(); i++)
        if(!record_processed_[i])
            return i;
    return size_t(-1);
}

std::unordered_set<size_t> Cdr3GraphComputer::GetRecordsThatShareCdr3sWith(size_t record_ind) {
    VERIFY_MSG(record_ind < pairing_data_cdr3s_.size(), "Index " << record_ind << " exceeds # complete records");
    std::unordered_set<size_t> indices;
    auto record_cdr3 = pairing_data_cdr3s_[record_ind];
    //INFO("#HCDR3s: " << record_cdr3.hcdr3s.size() <<
    //             ", # KCDR3s: " << record_cdr3.kcdr3s.size() <<
    //             ", # LCDR3s: " << record_cdr3.lcdr3s.size());
    for(auto it = record_cdr3.hcdr3s.begin(); it != record_cdr3.hcdr3s.end(); it++)
        indices.insert(hcdr3_map_[*it].begin(), hcdr3_map_[*it].end());
    for(auto it = record_cdr3.kcdr3s.begin(); it != record_cdr3.kcdr3s.end(); it++)
        indices.insert(kcdr3_map_[*it].begin(), kcdr3_map_[*it].end());
    for(auto it = record_cdr3.lcdr3s.begin(); it != record_cdr3.lcdr3s.end(); it++)
        indices.insert(lcdr3_map_[*it].begin(), lcdr3_map_[*it].end());
    return indices;
}

void Cdr3GraphComputer::ComputeCdr3GraphFromStartRecord(size_t start_record_ind) {
    std::queue<size_t> complete_record_queue;
    complete_record_queue.push(start_record_ind);
    std::unordered_set<size_t> related_records;
    //INFO("Index of start record: " << start_record_ind);
    while(!complete_record_queue.empty()) {
        size_t current_record = complete_record_queue.front();
        complete_record_queue.pop();
        related_records.insert(current_record);
        record_processed_[current_record] = true;
        std::unordered_set<size_t> indices = GetRecordsThatShareCdr3sWith(current_record);
        //INFO("# indices: " << indices.size());
        for (auto ind = indices.begin(); ind != indices.end(); ind++)
            if(!record_processed_[*ind]) {
                complete_record_queue.push(*ind);
                record_processed_[*ind] = true;
            }
    }
    //INFO("# related records: " << related_records.size());
    related_records_.push_back(related_records);
}

void Cdr3GraphComputer::ComputeCdr3Graphs() {
    InitializeMaps();
    size_t num_related_groups = 0;
    size_t start_record = GetStartRecordIndex();
    while(start_record != size_t(-1)) {
        ComputeCdr3GraphFromStartRecord(start_record);
        start_record = GetStartRecordIndex();
//        if(num_related_groups > 100)
//            break;
        num_related_groups++;
    }
}

void Cdr3GraphComputer::Compute() {
    ComputeHCdr3s();
    ComputeKCdr3s();
    ComputeLCdr3s();
    ComputeCdr3Graphs();
}

std::pair<size_t, size_t> Cdr3GraphComputer::ComputeIsotypesNumberForRelatedGroup(size_t index) {
    assert(index < related_records_.size());
    size_t num_lc_isotypes = 0;
    size_t num_kc_isotypes = 0;
    std::set<std::string> hc_isotypes;
    auto records = related_records_[index];
    for(auto it  = records.begin(); it != records.end(); it++) {
        auto cur_hc_isotypes = pairing_data_storage_[complete_indices_[*it]]->HcIsotypes();
        for(auto it2 = cur_hc_isotypes.begin(); it2 != cur_hc_isotypes.end(); it2++)
            hc_isotypes.insert(it2->str());
        if(pairing_data_storage_[complete_indices_[*it]]->KappaChainCount() > 0)
            num_kc_isotypes = 1;
        if(pairing_data_storage_[complete_indices_[*it]]->LambdaChainCount() > 0)
            num_lc_isotypes = 1;
    }
    return std::make_pair(hc_isotypes.size(), num_kc_isotypes + num_lc_isotypes);
}

std::string Cdr3GraphComputer::GetRelatedGroupDirName(size_t index, std::pair<size_t, size_t> num_isotypes) {
    assert(index < related_records_.size());
    std::stringstream ss;
    ss << "related_group_" << index << "_size_" << related_records_[index].size() <<
            "_hc_isotypes_" << num_isotypes.first << "_lc_isotypes_" << num_isotypes.second;
    return ss.str();
}

std::pair<size_t, size_t> Cdr3GraphComputer::ComputeNumberLightCdr3InRelatedGroup(size_t index) {
    assert(index < related_records_.size());
    std::unordered_set<Annotation, AnnotationHasher> kcdr3s;
    std::unordered_set<Annotation, AnnotationHasher> lcdr3s;
    for(auto it = related_records_[index].begin(); it != related_records_[index].end(); it++) {
        kcdr3s.insert(pairing_data_cdr3s_[*it].kcdr3s.begin(), pairing_data_cdr3s_[*it].kcdr3s.end());
        lcdr3s.insert(pairing_data_cdr3s_[*it].lcdr3s.begin(), pairing_data_cdr3s_[*it].lcdr3s.end());
    }
    return std::make_pair(kcdr3s.size(), lcdr3s.size());
}

bool Cdr3GraphComputer::OutputRelatedGroup(size_t index) {
    assert(index < related_records_.size());
    auto related_group = related_records_[index];
    if(related_group.size() > 30 or related_group.size() < 5)
        return false;
    auto num_isotypes = ComputeIsotypesNumberForRelatedGroup(index);
    //if(num_isotypes.second != 2)
    //    return false;
    return num_isotypes.second == 2;
    //auto num_light_cdr3 = ComputeNumberLightCdr3InRelatedGroup(index);
    //return num_light_cdr3.first > 1 or num_light_cdr3.second > 1;
}

void Cdr3GraphComputer::OutputRelatedGroupHcs(size_t index, std::string output_dir) {
    assert(index < related_records_.size());
    std::ofstream out(path::append_path(output_dir, "heavy_chains.fasta"));
    for(auto it = related_records_[index].begin(); it != related_records_[index].end(); it++) {
        size_t old_index = complete_indices_[*it];
        auto hc_isotypes = pairing_data_storage_[old_index]->HcIsotypes();
        for(auto hci = hc_isotypes.begin(); hci != hc_isotypes.end(); hci++) {
            auto seqs = pairing_data_storage_[old_index]->GetSequencesByIsotype(*hci);
            for(auto seq = seqs->cbegin(); seq != seqs->cend(); seq++) {
                out << ">" << GetHeaderForUmiSequence(pairing_data_storage_[old_index]->Db(), *seq) << std::endl;
                out << dna5string_to_stdstring(seq->sequence) << std::endl;
            }
        }
    }
}

void Cdr3GraphComputer::OutputRelatedGroupKcs(size_t index, std::string output_dir) {
    assert(index < related_records_.size());
    std::ofstream out(path::append_path(output_dir, "kappa_chains.fasta"));
    for(auto it = related_records_[index].begin(); it != related_records_[index].end(); it++) {
        size_t old_index = complete_indices_[*it];
        auto seqs = pairing_data_storage_[old_index]->GetSequencesByIsotype(IgIsotypeHelper::GetKappaIsotype());
        for(auto seq = seqs->cbegin(); seq != seqs->cend(); seq++) {
            out << ">" << GetHeaderForUmiSequence(pairing_data_storage_[old_index]->Db(), *seq) << std::endl;
            out << dna5string_to_stdstring(seq->sequence) << std::endl;
        }
    }
}

void Cdr3GraphComputer::OutputRelatedGroupLcs(size_t index, std::string output_dir) {
    assert(index < related_records_.size());
    std::ofstream out(path::append_path(output_dir, "lambda_chains.fasta"));
    for(auto it = related_records_[index].begin(); it != related_records_[index].end(); it++) {
        size_t old_index = complete_indices_[*it];
        auto seqs = pairing_data_storage_[old_index]->GetSequencesByIsotype(IgIsotypeHelper::GetLambdaIsotype());
        for(auto seq = seqs->cbegin(); seq != seqs->cend(); seq++) {
            out << ">" << GetHeaderForUmiSequence(pairing_data_storage_[old_index]->Db(), *seq) << std::endl;
            out << dna5string_to_stdstring(seq->sequence) << std::endl;
        }
    }
}

void Cdr3GraphComputer::OutputGraph(size_t index, std::string output_dir) {
    assert(index < related_records_.size());
    std::stringstream ss;
    ss << "graph_" << index << ".txt";
    std::ofstream out(path::append_path(output_dir, ss.str()));
    using std::endl;
    size_t num_db = 0;
    std::unordered_set<Annotation, AnnotationHasher> hcdr3s;
    std::unordered_set<Annotation, AnnotationHasher> kcdr3s;
    std::unordered_set<Annotation, AnnotationHasher> lcdr3s;
    // output droplet barcodes
    out << "Vertices:" << endl;
    for(auto it = related_records_[index].begin(); it != related_records_[index].end(); it++) {
        size_t old_index = complete_indices_[*it];
        out << num_db << "\t" <<
                pairing_data_storage_[old_index]->Db() << "\tgrey\trectangle" << endl;
        num_db++;
        hcdr3s.insert(pairing_data_cdr3s_[*it].hcdr3s.begin(), pairing_data_cdr3s_[*it].hcdr3s.end());
        kcdr3s.insert(pairing_data_cdr3s_[*it].kcdr3s.begin(), pairing_data_cdr3s_[*it].kcdr3s.end());
        lcdr3s.insert(pairing_data_cdr3s_[*it].lcdr3s.begin(), pairing_data_cdr3s_[*it].lcdr3s.end());
    }
    // output heavy chain cdr3s
    size_t num_hcdr3 = 0;
    std::unordered_map<Annotation, size_t, AnnotationHasher> hcdr3_ind;
    for(auto it = hcdr3s.begin(); it != hcdr3s.end(); it++) {
        out << num_db + num_hcdr3 << "\tH" << num_hcdr3 << ":" << it->v << "_" << it->j <<
                "\tyellow\tcircle" << endl;
        hcdr3_ind[*it] = num_db + num_hcdr3;
        num_hcdr3++;
    }
    // output kappa cdr3s
    size_t num_kcdr3 = 0;
    std::unordered_map<Annotation, size_t, AnnotationHasher> kcdr3_ind;
    for(auto it = kcdr3s.begin(); it != kcdr3s.end(); it++) {
        out << num_db + num_hcdr3 + num_kcdr3 << "\tK" << num_kcdr3 << ":" << it->v << "_" << it->j <<
                "\tred\tcircle" << endl;
        kcdr3_ind[*it] = num_db + num_hcdr3 + num_kcdr3;
        num_kcdr3++;
    }
    // output lambda cdr3s
    size_t num_lcdr3 = 0;
    std::unordered_map<Annotation, size_t, AnnotationHasher> lcdr3_ind;
    for(auto it = lcdr3s.begin(); it != lcdr3s.end(); it++) {
        out << num_db + num_hcdr3 + num_kcdr3 + num_lcdr3 << "\tL" << num_lcdr3 << ":" << it->v << "_" << it->j <<
                "\tgreen\tcircle" << endl;
        lcdr3_ind[*it] = num_db + num_hcdr3 + num_kcdr3 + num_lcdr3;
        num_lcdr3++;
    }
    // output edges
    out << "Edges:" << endl;
    size_t vindex = 0;
    for(auto ind = related_records_[index].begin(); ind != related_records_[index].end(); ind++) {
        for(auto it = pairing_data_cdr3s_[*ind].hcdr3s.begin(); it != pairing_data_cdr3s_[*ind].hcdr3s.end(); it++)
            out << vindex << "\t" << hcdr3_ind[*it] << endl;
        for(auto it = pairing_data_cdr3s_[*ind].kcdr3s.begin(); it != pairing_data_cdr3s_[*ind].kcdr3s.end(); it++)
            out << vindex << "\t" << kcdr3_ind[*it] << endl;
        for(auto it = pairing_data_cdr3s_[*ind].lcdr3s.begin(); it != pairing_data_cdr3s_[*ind].lcdr3s.end(); it++)
            out << vindex << "\t" << lcdr3_ind[*it] << endl;
        vindex++;
    }
    out.close();
}

void Cdr3GraphComputer::OutputCdr3Graphs() {
    if(!io_config_.output.output_related_groups)
        return;
    size_t num_reported = 0;
    for(size_t i = 0; i < related_records_.size(); i++) {
        if(!OutputRelatedGroup(i))
            continue;
        auto num_isotypes = ComputeIsotypesNumberForRelatedGroup(i);
        auto related_group = related_records_[i];
        INFO("Related group " << i << " (size: " << related_group.size() << "):");
        INFO("# HC isotypes: " << num_isotypes.first << ", # LC isotypes: " << num_isotypes.second);
        std::string output_dir = path::append_path(io_config_.output.related_groups_dir,
                                                   GetRelatedGroupDirName(i, num_isotypes));
        path::make_dir(output_dir);
        OutputRelatedGroupHcs(i, output_dir);
        OutputRelatedGroupKcs(i, output_dir);
        OutputRelatedGroupLcs(i, output_dir);
        OutputGraph(i, output_dir);
        INFO("Sequences were written to " << output_dir);
        //for(auto it2 = related_group.begin(); it2 != related_group.end(); it2++)
        //    std::cout << *(pairing_data_storage_[complete_indices_[*it2]]) << std::endl;
        INFO("--------------------");
        num_reported++;
    }
    INFO("# reported groups: " << num_reported);
}
