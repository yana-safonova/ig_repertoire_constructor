//
// Created by Andrew Bzikadze on 5/18/16.
//

#include <tuple>
#include <string>

#include <seqan/seq_io.h>

#include "alignment_reader.hpp"

namespace shm_kmer_matrix_estimator {

AlignmentReader::AlignmentReader(const std::string &alignments_filename,
                                 const std::string &cdr_details_filename,
                                 const shm_kmer_matrix_estimator_config::alignment_checker_params &alignment_checker_params,
                                 const shm_kmer_matrix_estimator_config::alignment_cropper_params &alignment_cropper_params) :
    alignments_filename_(alignments_filename),
    cdr_details_filename_(cdr_details_filename)
{
    using AlignmentCheckerMethod = shm_kmer_matrix_estimator_config::alignment_checker_params::AlignmentCheckerMethod;
    if (alignment_checker_params.alignment_checker_method == AlignmentCheckerMethod::NoGaps) {
        alignment_checker_ptr_ = std::unique_ptr<NoGapsAlignmentChecker>
            (new NoGapsAlignmentChecker(alignment_checker_params));
    }

    using AlignmentCropperMethod = shm_kmer_matrix_estimator_config::alignment_cropper_params::AlignmentCropperMethod;
    if (alignment_cropper_params.alignment_cropper_method == AlignmentCropperMethod::UptoLastReliableKMer) {
        alignment_cropper_ptr_ = std::unique_ptr<UptoLastReliableKmerAlignmentCropper>
            (new UptoLastReliableKmerAlignmentCropper(alignment_cropper_params.rkmp));
    }
}

VectorEvolutionaryEdgeAlignments AlignmentReader::read_alignments() const {
    // cdr_details
    std::ifstream cdr_details(cdr_details_filename_);
    VERIFY_MSG(cdr_details, "File is not open");
    std::string cdr_details_line;
    std::getline(cdr_details, cdr_details_line); // skip header

    // v_alignments
    std::vector<seqan::CharString> names;
    std::vector<seqan::CharString> reads;
    seqan::SeqFileIn seqFileIn(alignments_filename_.c_str());
    seqan::readRecords(names, reads, seqFileIn);
    VectorEvolutionaryEdgeAlignments alignments;

    auto ReadIterator = reads.cbegin();
    auto NamesIterator = names.cbegin();

    while (ReadIterator != reads.cend()) {
        // Reading v_alignment
        std::string read_seq = std::string(seqan::toCString(*ReadIterator));
        ++ReadIterator;
        ++NamesIterator;
        VERIFY(ReadIterator != reads.cend());
        VERIFY(NamesIterator != names.cend());
        std::string germline_seq = std::string(seqan::toCString(*ReadIterator));
        std::string gene_id = std::string(seqan::toCString(*NamesIterator));
        ++ReadIterator;
        ++NamesIterator;

        // Reading cdr info
        std::getline(cdr_details, cdr_details_line);
        std::istringstream cdr_details_stream(cdr_details_line);
        std::string temp;
        bool has_stop_codon, in_frame, productive;
        size_t cdr1_start, cdr1_end, cdr2_start, cdr2_end;
        for (size_t i = 0; i < has_stop_codon_pos; ++i) {
            cdr_details_stream >> temp;
        }
        cdr_details_stream >> has_stop_codon >> in_frame >> productive >> temp >>
                 cdr1_start >> cdr1_end >> temp >> cdr2_start >> cdr2_end;
        // if (not (cdr_details_stream))
        //     continue;

        --cdr1_start, --cdr1_end, --cdr2_start, --cdr2_end;
        EvolutionaryEdgeAlignment alignment(std::move(germline_seq), std::move(read_seq), gene_id,
                                            has_stop_codon, in_frame, productive,
                                            cdr1_start, cdr1_end, cdr2_start, cdr2_end);

        if (alignment_checker_ptr_->check(alignment)) {
            alignment_cropper_ptr_->crop(alignment);
            alignments.emplace_back(std::move(alignment));
        }
    }
    return alignments;
}

} // End namespace shm_kmer_matrix_estimator
