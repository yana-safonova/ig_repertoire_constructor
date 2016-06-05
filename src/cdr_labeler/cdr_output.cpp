#include "cdr_output.hpp"

#include <seqan/seq_io.h>
#include <seqan/stream.h>

namespace cdr_labeler {
    std::ostream& CDRLabelingWriter::OutputCloneRegion(std::ostream& out,
                                                       const annotation_utils::CDRAnnotatedClone &clone,
                                                       annotation_utils::StructuralRegion region) const {
        if(clone.RegionIsEmpty(region)) {
            out << "-\t-\t-";
            return out;
        }
        auto region_range = clone.GetRangeByRegion(region);
        out << clone.GetRegionString(region) << "\t" << region_range.start_pos + 1 << "\t" << region_range.end_pos + 1;
        return out;
    }

    void CDRLabelingWriter::OutputCDRDetails() const {
        std::ofstream out(output_config_.cdr_details);
        out << "Read_name\tV_hit\tJ_hit\tCDR1_nucls\tCDR1_start\tCDR1_end\tCDR2_nucls\tCDR2_start\tCDR2_end\tCDR3_nucls\tCDR3_start\tCDR3_end" << std::endl;
        for(auto it = clone_set_.cbegin(); it != clone_set_.cend(); it++) {
            annotation_utils::CDRAnnotatedClone cdr_clone = *it;
            auto vj_hit = alignment_info_.GetVJHitsByRead(cdr_clone.Read());
            out << cdr_clone.Read().name << "\t" << vj_hit.GetVHitByIndex(0).ImmuneGene().name() << "\t" <<
                    vj_hit.GetJHitByIndex(0).ImmuneGene().name() << "\t";
            OutputCloneRegion(out, cdr_clone, annotation_utils::StructuralRegion::CDR1);
            out << "\t";
            OutputCloneRegion(out, cdr_clone, annotation_utils::StructuralRegion::CDR2);
            out << "\t";
            OutputCloneRegion(out, cdr_clone, annotation_utils::StructuralRegion::CDR3);
            out << std::endl;
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