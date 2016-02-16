#include "logger/logger.hpp"

#include "raw_pairing_data.hpp"
#include "pairing_fastq_utils.hpp"

void RawPairingData::CheckDropletBarcodeConsistency(DropletBarcode db) {
    assert(droplet_barcode_ == db);
}

void RawPairingData::UpdateHcSequences(IsotypeUmiSequence hc_sequence) {
    if(hc_isotype_map_.find(hc_sequence.isotype) == hc_isotype_map_.end())
        hc_isotype_map_[hc_sequence.isotype] = IsotypeUmiSequencesPtr(new IsotypeUmiSequences(hc_sequence.isotype));
    hc_isotype_map_[hc_sequence.isotype]->Update(hc_sequence);
}

void RawPairingData::Update(DropletBarcode db, std::string header, std::string sequence) {
    CheckDropletBarcodeConsistency(db);
    IsotypeUmiSequence umi_sequence(PairingFastqUtils::ExtractIsotypeFromHeader(header),
                                    PairingFastqUtils::ExtractUmiFromHeader(header),
                                    PairingFastqUtils::ExtractSizeFromHeader(header),
                                    PairingFastqUtils::ConvertToDnaString(sequence));
    if(umi_sequence.isotype.IsHeavyChain()) {
        UpdateHcSequences(umi_sequence);
    }
    else if(umi_sequence.isotype.IsIgK()) {
        kappa_sequences_->Update(umi_sequence);
    }
    else if(umi_sequence.isotype.IsIgL()) {
        lambda_sequences_->Update(umi_sequence);
    }
}

std::vector<IgIsotype> RawPairingData::HcIsotypes() {
    if(hc_isotypes_.size() == 0)
        for(auto it = hc_isotype_map_.begin(); it != hc_isotype_map_.end(); it++)
            hc_isotypes_.push_back(it->first);
    return hc_isotypes_;
}

bool RawPairingData::ContainsHcIsotype(IgIsotype isotype) const {
    for(auto it = hc_isotypes_.cbegin(); it != hc_isotypes_.cend() ;it++)
        if(isotype == *it)
            return true;
    return false;
}

size_t RawPairingData::HcRecordsCount(IgIsotype hc_isotype) const {
    if(hc_isotype_map_.find(hc_isotype) == hc_isotype_map_.end())
        return 0;
    return hc_isotype_map_.at(hc_isotype)->size();
}

IsotypeUmiSequencesPtr RawPairingData::GetSequencesByIsotype(IgIsotype isotype) const {
    if(isotype.IsIgK())
        return kappa_sequences_;
    if(isotype.IsIgL())
        return lambda_sequences_;
    assert(isotype.IsHeavyChain());
    if(ContainsHcIsotype(isotype))
        return hc_isotype_map_.at(isotype);
    return IsotypeUmiSequencesPtr(new IsotypeUmiSequences());
}

size_t RawPairingData::TotalNumberHcs() const {
    size_t total_num_hcs = 0;
    for(auto it = hc_isotype_map_.cbegin(); it != hc_isotype_map_.cend(); it++)
        total_num_hcs += it->second->size();
    return total_num_hcs;
}