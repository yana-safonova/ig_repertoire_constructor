#pragma once

#include <abpair_config.hpp>
#include "raw_pairing_data_storage.hpp"

class RawPairingDataStatsCalculator {
    // input
    const abpair_config::io_config::output_config &output_;
    const RawPairingDataStorage &raw_pairing_data_storage_;

    // statistics to be computed
    size_t num_complete_records_;
    size_t num_complete_non_ambiguous_records_;
    size_t num_complete_non_ambiguous_hc_records_;
    size_t num_hc_ambiguous_records_;
    size_t num_both_lc_records_;
    size_t num_hc_lc_ambiguous_records_;
//    size_t num_kappa_records_;
//    size_t num_lambda_records_;
//    size_t num_igg_records_;
//    size_t num_igm_records_;
//    size_t num_igd_records_;
//    size_t num_iga_records_;
//    size_t num_ige_records_;

    void ComputeStats();

    std::string GetFilenameForRawRecord(RawPairingDataPtr pairing_record, std::string isotype_str);

    std::string GetHeaderForUmiSequence(const IsotypeUmiSequence &umi_sequence);

    std::string GetBarcodeDir(RawPairingDataPtr pairing_record) const;

    void OutputBarcodesByIsotype(IsotypeUmiSequencesPtr umi_sequences, std::string barcode_dir, std::string db);

    std::string GetOutputFnameForIsotypeBarcodes(IgIsotype isotype, size_t size);

    bool BarcodeShouldBeReported(RawPairingDataPtr pairing_record);

public:
    RawPairingDataStatsCalculator(const abpair_config::io_config::output_config &output,
                                  const RawPairingDataStorage &raw_pairing_data_storage) :
            output_(output),
            raw_pairing_data_storage_(raw_pairing_data_storage),
            num_complete_records_(),
            num_complete_non_ambiguous_records_(),
            num_complete_non_ambiguous_hc_records_(),
            num_hc_ambiguous_records_(),
            num_both_lc_records_(),
            num_hc_lc_ambiguous_records_()
//            num_kappa_records_(),
//            num_lambda_records_(),
//            num_igg_records_(),
//            num_igm_records_(),
//            num_igd_records_(),
//            num_iga_records_(),
//            num_ige_records_()
    {
        ComputeStats();
    }

    size_t NumberCompleteRecords() const { return num_complete_records_; }

    size_t NumberCompleteNonAmbiguousRecords() const { return num_complete_non_ambiguous_records_; }

    size_t NumberCompleteNonAmbiguousHcRecords() const { return num_complete_non_ambiguous_hc_records_; }

    size_t NumberHcAmbiguousRecords() const { return num_hc_ambiguous_records_; }

    size_t NumberRecordsWithBothLcs() const { return num_both_lc_records_; }

    size_t NumberHcLcAmbiguous() const { return num_hc_lc_ambiguous_records_; }

    void OutputHcAmbiguousRecords();

    void OutputMolecularBarcodes();

    void OutputDemultiplexedData();
};