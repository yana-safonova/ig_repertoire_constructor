#include "logger/logger.hpp"

#include "raw_pairing_data.hpp"
#include "pairing_fastq_utils.hpp"

std::string dna5string_to_stdstring(seqan::Dna5String dna5string) {
    seqan::String<char> char_string;
    seqan::assign(char_string, dna5string);
    return std::string(seqan::toCString(char_string));
}

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

std::ostream& operator<<(std::ostream& out, RawPairingData &raw_pairing_data) {
    out << "Droplet barcode: " << raw_pairing_data.Db() << std::endl;
    if(!raw_pairing_data.HcIsMissed()) {
        out << "# HC isotypes: " << raw_pairing_data.HcIsotypes().size() << std::endl;
        auto hc_isotypes = raw_pairing_data.HcIsotypes();
        for (auto hc = hc_isotypes.begin(); hc != hc_isotypes.end(); hc++) {
            auto hc_seqs = raw_pairing_data.GetSequencesByIsotype(*hc);
            out << "HC isotype: " << hc->str() << ", # sequences: " << hc_seqs->size() << std::endl;
            for (auto hc_seq = hc_seqs->cbegin(); hc_seq != hc_seqs->cend(); hc_seq++)
                out << dna5string_to_stdstring(hc_seq->sequence) << std::endl;
        }
    }
    if(raw_pairing_data.KappaChainCount() != 0) {
        auto kappa_seqs = raw_pairing_data.GetSequencesByIsotype(IgIsotypeHelper::GetKappaIsotype());
        out << "KC isotype, # sequences: " << kappa_seqs->size() << std::endl;
        for(auto it = kappa_seqs->cbegin(); it != kappa_seqs->cend(); it++)
            out << dna5string_to_stdstring(it->sequence) << std::endl;
    }
    if(raw_pairing_data.LambdaChainCount() != 0) {
        auto lambda_seqs = raw_pairing_data.GetSequencesByIsotype(IgIsotypeHelper::GetLambdaIsotype());
        out << "LC isotype, # sequences: " << lambda_seqs->size() << std::endl;
        for(auto it = lambda_seqs->cbegin(); it != lambda_seqs->cend(); it++)
            out << dna5string_to_stdstring(it->sequence) << std::endl;
    }
    return out;
}
