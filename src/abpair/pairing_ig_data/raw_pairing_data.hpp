#pragma once

#include "pairing_primitives/ig_isotype.hpp"
#include "pairing_primitives/droplet_barcode.hpp"
#include "pairing_primitives/umi_isotype_sequences.hpp"

std::string dna5string_to_stdstring(seqan::Dna5String dna5string);

class RawPairingData {
    DropletBarcode droplet_barcode_;

    std::map<IgIsotype, IsotypeUmiSequencesPtr> hc_isotype_map_;
    std::vector<IgIsotype> hc_isotypes_;

    IsotypeUmiSequencesPtr kappa_sequences_;
    IsotypeUmiSequencesPtr lambda_sequences_;

    void CheckDropletBarcodeConsistency(DropletBarcode db);

    void UpdateHcSequences(IsotypeUmiSequence hc_sequence);

public:
    RawPairingData(DropletBarcode db) :
            droplet_barcode_(db),
            kappa_sequences_(new IsotypeUmiSequences(IgIsotypeHelper::GetKappaIsotype())),
            lambda_sequences_(new IsotypeUmiSequences(IgIsotypeHelper::GetLambdaIsotype())) { }

    void Update(DropletBarcode db, std::string header, std::string sequence,
                std::string v_gene, std::string j_gene);

    DropletBarcode Db() const { return droplet_barcode_; }

    size_t LambdaChainCount() const { return lambda_sequences_->size(); }

    size_t KappaChainCount() const { return kappa_sequences_->size(); }

    size_t HcIsotypeNumber() const { return hc_isotype_map_.size(); }

    std::vector<IgIsotype> HcIsotypes();

    bool ContainsHcIsotype(IgIsotype isotype) const;

    size_t HcRecordsCount(IgIsotype hc_isotype) const;

    IsotypeUmiSequencesPtr GetSequencesByIsotype(IgIsotype isotype) const;

    bool HcIsAmbiguous() const { return HcIsotypeNumber() > 1; }

    bool HcIsNonAmbiguous() const { return HcIsotypeNumber() == 1; }

    bool LcIsAmbiguous() const { return LambdaChainCount() > 0 and KappaChainCount() > 0; }

    bool LcIsNonAmbiguous() const { return (KappaChainCount() > 0 and LambdaChainCount() == 0) or
                (KappaChainCount() == 0 and LambdaChainCount() > 0); }

    bool HcIsMissed() const { return HcIsotypeNumber() == 0; }

    bool LcIsMissed() const { return KappaChainCount() == 0 and LambdaChainCount() == 0; }

    bool Complete() const { return !HcIsMissed() and !LcIsMissed(); }

    bool CompleteAndNonHcAmbiguous() const { return Complete() and !HcIsAmbiguous(); }

    size_t TotalNumberHcs() const;
};

std::ostream& operator<<(std::ostream&, RawPairingData &raw_pairing_data);

typedef std::shared_ptr<RawPairingData> RawPairingDataPtr;