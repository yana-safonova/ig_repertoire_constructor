#pragma once

#include "ig_isotype.hpp"
#include "seqan/sequence.h"

struct UmiIsotypeSequence {
    IgIsotype isotype;
    std::string umi;
    size_t size;
    seqan::Dna5String sequence;

    UmiIsotypeSequence() :
            isotype(),
            umi(),
            size(),
            sequence() { }

    UmiIsotypeSequence(IgIsotype new_isotype,
                       std::string new_umi,
                       size_t new_size,
                       seqan::Dna5String new_sequence) :
            isotype(new_isotype),
            umi(new_umi),
            size(new_size),
            sequence(new_sequence) { }

};

//--------------------------------------------------------------

class IsotypeUmiSequences {
    IgIsotype isotype_;
    std::vector<UmiIsotypeSequence> sequences_;

    void CheckConsistency(const UmiIsotypeSequence &umi_sequence) const;

public:
    IsotypeUmiSequences() : isotype_() { }

    IsotypeUmiSequences(IgIsotype isotype) : isotype_(isotype) { }

    void Update(UmiIsotypeSequence umi_sequence);

    size_t size() const { return sequences_.size(); }

    typedef std::vector<UmiIsotypeSequence>::const_iterator umi_seq_citerator;

    umi_seq_citerator cbegin() const { return sequences_.cbegin(); }

    umi_seq_citerator cend() const { return sequences_.cend(); }

    IgIsotype Isotype() const { return isotype_; }
};

typedef std::shared_ptr<IsotypeUmiSequences> IsotypeUmiSequencesPtr;

//--------------------------------------------------------------

class RawPairingData {
    std::string droplet_barcode_;

    std::map<IgIsotype, IsotypeUmiSequencesPtr> hc_isotype_map_;
    std::vector<IgIsotype> hc_isotypes_;

    IsotypeUmiSequencesPtr kappa_sequences_;
    IsotypeUmiSequencesPtr lambda_sequences_;

    void CheckDropletBarcodeConsistency(std::string db);

    void UpdateHcSequences(UmiIsotypeSequence hc_sequence);

public:
    RawPairingData(std::string db) :
            droplet_barcode_(db),
            kappa_sequences_(new IsotypeUmiSequences(IgIsotypeHelper::GetKappaIsotype())),
            lambda_sequences_(new IsotypeUmiSequences(IgIsotypeHelper::GetLambdaIsotype())) { }

    void Update(std::string header, std::string sequence);

    std::string Db() const { return droplet_barcode_; }

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

typedef std::shared_ptr<RawPairingData> RawPairingDataPtr;