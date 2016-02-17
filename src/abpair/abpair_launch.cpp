#include "logger/logger.hpp"
#include "verify.hpp"
#include "abpair_launch.hpp"
#include "pairing_ig_data/raw_pairing_data_storage.hpp"
#include "pairing_ig_data/raw_pairing_data_stats_calculator.hpp"
#include "cdr3_computation/simple_cdr3_calculator.hpp"

std::vector<std::string> AbPairLauncher::ReadInputFnames(std::string input_sequences) {
    std::ifstream ifhandler(input_sequences);
    VERIFY_MSG(ifhandler.good(), "File " << input_sequences << " was not found");
    std::vector<std::string> fnames;
    while(!ifhandler.eof()) {
        std::string tmp;
        std::getline(ifhandler, tmp);
        if(tmp != "")
            fnames.push_back(tmp);
    }
    return fnames;
}

void AbPairLauncher::Run(const abpair_config::io_config &io) {
    INFO("==== AbPair starts");
    // reading input sequences
    std::vector<std::string> fastq_fnames = ReadInputFnames(io.input.input_sequences);
    INFO("Reading input sequences");
    INFO(fastq_fnames.size() << " FASTQ files will be extracted from " << io.input.input_sequences);
    RawPairingDataStorage raw_pairing_storage;
    for(auto it = fastq_fnames.begin(); it != fastq_fnames.end(); it++)
        raw_pairing_storage.Update(*it);
    INFO(raw_pairing_storage.size() << " pairing records were extracted from input reads");
    //for(auto it = raw_pairing_storage.cbegin(); it != raw_pairing_storage.cend(); it++) {
    //    std::cout << (*it)->Db() << ". #HC isotypes: " << (*it)->HcIsotypeNumber() << ", #Ks: " <<
    //            (*it)->KappaChainCount() << ", #Ls: " << (*it)->LambdaChainCount() << std::endl;
    //}
    RawPairingDataStatsCalculator stats_calculator(io.output, raw_pairing_storage);
    INFO("# complete records: " << stats_calculator.NumberCompleteRecords());
    INFO("# complete and non-ambiguous records: " << stats_calculator.NumberCompleteNonAmbiguousRecords());
    INFO("# complete and non-ambiguous HC records: " << stats_calculator.NumberCompleteNonAmbiguousHcRecords());
    INFO("# HC ambiguous records: " << stats_calculator.NumberHcAmbiguousRecords());
    INFO("# records with KC & LC: " << stats_calculator.NumberRecordsWithBothLcs());
    INFO("# records with ambiguous HCs and LCs: " << stats_calculator.NumberHcLcAmbiguous());
    stats_calculator.OutputHcAmbiguousRecords();
    stats_calculator.OutputMolecularBarcodes();
    stats_calculator.OutputDemultiplexedData();
    INFO("CDR3 computation starts");
    for(auto it = raw_pairing_storage.cbegin(); it != raw_pairing_storage.cend(); it++) {
        if(!(*it)->Complete())
            continue;
        INFO("IGH CDR3");
        auto hc_isotypes = (*it)->HcIsotypes();
        for(auto hc = hc_isotypes.begin(); hc != hc_isotypes.end(); hc++) {
            auto hc_seqs = (*it)->GetSequencesByIsotype(*hc);
            for(auto hc_seq = hc_seqs->cbegin(); hc_seq != hc_seqs->cend(); hc_seq++)
                auto pos = SimpleCdr3Calculator().FindCdr3Positions((*hc_seq).sequence);
        }
        if((*it)->KappaChainCount() != 0) {
            INFO("IGK CDR3");
            auto kappa_seq = (*it)->GetSequencesByIsotype(IgIsotypeHelper::GetKappaIsotype());
            for(auto kit = kappa_seq->cbegin(); kit != kappa_seq->cend(); kit++)
                auto pos = SimpleCdr3Calculator().FindCdr3Positions((*kit).sequence);
        }
        if((*it)->LambdaChainCount() != 0) {
            INFO("IGL CDR3");
            auto lambda_seq = (*it)->GetSequencesByIsotype(IgIsotypeHelper::GetLambdaIsotype());
            for(auto lit = lambda_seq->cbegin(); lit != lambda_seq->cend(); lit++)
                auto pos = SimpleCdr3Calculator().FindCdr3Positions((*lit).sequence);
        }
        break;
    }
    INFO("==== AbPair ends");
}