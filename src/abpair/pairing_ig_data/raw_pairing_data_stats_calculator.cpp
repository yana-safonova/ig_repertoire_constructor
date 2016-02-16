#include "logger/logger.hpp"
#include "raw_pairing_data_stats_calculator.hpp"
#include <sstream>
#include <seqan/seq_io.h>

void RawPairingDataStatsCalculator::ComputeStats() {
    for(auto it = raw_pairing_data_storage_.cbegin(); it != raw_pairing_data_storage_.cend(); it++) {
        if((*it)->Complete())
            num_complete_records_++;
        if((*it)->CompleteAndNonHcAmbiguous())
            num_complete_non_ambiguous_hc_records_++;
        if((*it)->CompleteAndNonHcAmbiguous() and (*it)->LcIsNonAmbiguous())
            num_complete_non_ambiguous_records_++;
        if((*it)->HcIsAmbiguous())
            num_hc_ambiguous_records_++;
        if((*it)->LcIsAmbiguous())
            num_both_lc_records_++;
        if((*it)->HcIsAmbiguous() and (*it)->LcIsAmbiguous())
            num_hc_lc_ambiguous_records_++;
    }
}

std::string RawPairingDataStatsCalculator::GetFilenameForRawRecord(RawPairingDataPtr pairing_record) {
    std::stringstream ss;
    ss << "db_" << pairing_record->Db() << "_isotypes_" << pairing_record->HcIsotypeNumber() <<
            "_seqs_" << pairing_record->TotalNumberHcs() << ".fasta";
    return path::append_path(output_.hc_ambiguous_dir, ss.str());
}

std::string RawPairingDataStatsCalculator::GetHeaderForUmiSequence(const IsotypeUmiSequence &umi_sequence) {
    return std::string("UMI:" + umi_sequence.umi + "|ISOTYPE:" + umi_sequence.isotype.str());
}

void RawPairingDataStatsCalculator::OutputHcAmbiguousRecords() {
    INFO("Ambiguous heavy chain sequences will be written to " + output_.hc_ambiguous_dir);
    size_t num_written_files = 0;
    for(auto it = raw_pairing_data_storage_.cbegin(); it != raw_pairing_data_storage_.cend(); it++) {
        if((*it)->HcIsNonAmbiguous() or (*it)->HcIsMissed())
            continue;
        auto hc_isotypes = (*it)->HcIsotypes();
        seqan::SeqFileOut out(GetFilenameForRawRecord(*it).c_str());
        num_written_files++;
        for(auto isotype = hc_isotypes.begin(); isotype != hc_isotypes.end(); isotype++) {
            auto umi_sequences = (*it)->GetSequencesByIsotype(*isotype);
            for(auto umi_seq = umi_sequences->cbegin(); umi_seq != umi_sequences->cend(); umi_seq++) {
                seqan::writeRecord(out,
                                   seqan::CharString(GetHeaderForUmiSequence(*umi_seq)),
                                   umi_seq->sequence);
            }
        }
    }
    INFO(num_written_files << " files were written to " << output_.hc_ambiguous_dir);
}

std::string RawPairingDataStatsCalculator::GetBarcodeDir(RawPairingDataPtr pairing_record) const {
    std::stringstream ss;
    ss << pairing_record->Db();
    return path::append_path(output_.barcode_dir, ss.str());
}

std::string RawPairingDataStatsCalculator::GetOutputFnameForIsotypeBarcodes(IgIsotype isotype,
                                                                            size_t size) {
    std::stringstream ss;
    ss << "isotype_" << isotype.str() << "_size_" << size << ".fasta";
    return ss.str();
}

void RawPairingDataStatsCalculator::OutputBarcodesByIsotype(IsotypeUmiSequencesPtr umi_sequences,
                                                            std::string barcode_dir,
                                                            std::string db) {
    if(umi_sequences->size() == 0)
        return;
    std::string output_fname = path::append_path(barcode_dir,
                                                 GetOutputFnameForIsotypeBarcodes(umi_sequences->Isotype(),
                                                                                  umi_sequences->size()));
    std::ofstream output_fhandler(output_fname);
    for(auto it = umi_sequences->cbegin(); it != umi_sequences->cend(); it++) {
        output_fhandler << ">DB:" << db <<  "|MB:" << it->umi << "|ISOTYPE:" << it->isotype.str() <<
                "|COUNT:" << it->size << std::endl;
        output_fhandler << it->umi << std::endl;
    }
    output_fhandler.close();
}

bool RawPairingDataStatsCalculator::BarcodeShouldBeReported(RawPairingDataPtr pairing_record) {
    return true;
    auto hc_isotypes = pairing_record->HcIsotypes();
    for(auto it = hc_isotypes.begin(); it != hc_isotypes.end(); it++)
        if(pairing_record->GetSequencesByIsotype(*it)->size() < 2)
            return false;
    if(pairing_record->GetSequencesByIsotype(IgIsotypeHelper::GetKappaIsotype())->size() < 2)
        return false;
    return pairing_record->GetSequencesByIsotype(IgIsotypeHelper::GetLambdaIsotype())->size() > 1;
}

void RawPairingDataStatsCalculator::OutputMolecularBarcodes() {
    INFO("Molecular barcodes will be written to " + output_.barcode_dir);
    size_t num_reported = 0;
    for(auto it = raw_pairing_data_storage_.cbegin(); it != raw_pairing_data_storage_.cend(); it++) {
        if(!BarcodeShouldBeReported(*it))
            continue;
        std::string barcode_dir = GetBarcodeDir(*it);
        path::make_dir(barcode_dir);
        std::string db = (*it)->Db().StrId();
        OutputBarcodesByIsotype((*it)->GetSequencesByIsotype(IgIsotypeHelper::GetKappaIsotype()), barcode_dir, db);
        OutputBarcodesByIsotype((*it)->GetSequencesByIsotype(IgIsotypeHelper::GetLambdaIsotype()), barcode_dir, db);
        auto hc_isotypes = (*it)->HcIsotypes();
        for(auto ig_it = hc_isotypes.begin(); ig_it != hc_isotypes.end(); ig_it++)
            OutputBarcodesByIsotype((*it)->GetSequencesByIsotype(*ig_it), barcode_dir, db);
        num_reported++;
    }
    INFO(num_reported << " droplet barcodes were reported");
}