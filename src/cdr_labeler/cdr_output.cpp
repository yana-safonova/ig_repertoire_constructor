#include "cdr_output.hpp"

#include <seqan/seq_io.h>
#include <seqan/stream.h>

namespace cdr_labeler {
    void CDRLabelingWriter::OutputCDRDetails() const {
        std::ofstream out(output_config_.cdr_details);
        out << "Read_name\tV_hit\tJ_hit\tCDR1_nucls\tCDR1_start\tCDR1_end\tCDR2_nucls\tCDR2_start\tCDR2_end\tCDR3_nucls\tCDR3_start\tCDR3_end" << std::endl;
        for(auto it = clone_set_.cbegin(); it != clone_set_.cend(); it++) {
            annotation_utils::CDRAnnotatedClone cdr_clone = *it;
            auto cdr1_range = cdr_clone.GetRangeByRegion(annotation_utils::StructuralRegion::CDR1);
            auto cdr2_range = cdr_clone.GetRangeByRegion(annotation_utils::StructuralRegion::CDR2);
            auto cdr3_range = cdr_clone.GetRangeByRegion(annotation_utils::StructuralRegion::CDR3);
            auto vj_hit = alignment_info_.GetVJHitsByRead(cdr_clone.Read());
            out << cdr_clone.Read().name << "\t" << vj_hit.GetVHitByIndex(0).ImmuneGene().name() << "\t" <<
                    vj_hit.GetJHitByIndex(0).ImmuneGene().name() << "\t" << cdr_clone.CDR1() << "\t" <<
                    cdr1_range.start_pos + 1 << "\t" << cdr1_range.end_pos + 1 << "\t" <<
                    cdr_clone.CDR2() << "\t" << cdr2_range.start_pos + 1 << "\t" <<
                    cdr2_range.end_pos + 1 << "\t" << cdr_clone.CDR3() << "\t" <<
                    cdr3_range.start_pos + 1 << "\t" << cdr3_range.end_pos + 1 << std::endl;
        }
        out.close();
        INFO("CDR details were written to " << output_config_.cdr_details);
    }

    void CDRLabelingWriter::OutputCDR3Fasta() const {
        seqan::SeqFileOut out(output_config_.cdr3_fasta.c_str());
        std::vector<seqan::CharString> headers;
        std::vector<seqan::Dna5String> cdr3s;
        for(auto it = clone_set_.cbegin(); it != clone_set_.cend(); it++) {
            headers.push_back(seqan::CharString(it->Read().name.c_str()));
            cdr3s.push_back(it->CDR3());
        }
        seqan::writeRecords(out, headers, cdr3s);
        INFO("CDR3 sequences were written to " << output_config_.cdr3_fasta);
    }
}