#include "cdr_output.hpp"

#include "compressed_cdr_set.hpp"

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

    void CDRLabelingWriter::OutputRegionFasta(std::string output_fname,
                                              annotation_utils::StructuralRegion region) const {
        seqan::SeqFileOut out(output_fname.c_str());
        std::vector<seqan::CharString> headers;
        std::vector<seqan::Dna5String> regions;
        for(auto it = clone_set_.cbegin(); it != clone_set_.cend(); it++) {
            if(it->RegionIsEmpty(region))
                continue;
            headers.push_back(seqan::CharString(it->Read().name.c_str()));
            regions.push_back(it->GetRegionString(region));
        }
        seqan::writeRecords(out, headers, regions);
        INFO(region << " sequences were written to " << output_fname);
    }

    void CDRLabelingWriter::OutputCDR1Fasta() const {
        OutputRegionFasta(output_config_.cdr1_fasta, annotation_utils::StructuralRegion::CDR1);
    }

    void CDRLabelingWriter::OutputCDR2Fasta() const {
        OutputRegionFasta(output_config_.cdr2_fasta, annotation_utils::StructuralRegion::CDR2);
    }

    void CDRLabelingWriter::OutputCDR3Fasta() const {
        OutputRegionFasta(output_config_.cdr3_fasta, annotation_utils::StructuralRegion::CDR3);
    }

    seqan::CharString CDRLabelingWriter::GetCompressedRegionFname(annotation_utils::StructuralRegion region,
                                                                  CDRKey cdr_key, size_t abundance) const {
        std::stringstream ss;
        ss << region << ":" << cdr_key.id + 1 << "|V_hit:" << cdr_key.v_name << "|J_hit:" << cdr_key.j_name <<
                "|COUNT:" << abundance;
        return seqan::CharString(ss.str().c_str());
    }

    void CDRLabelingWriter::OutputCompressedCDR3Fasta() const {
        CompressedCDRSet compressed_cdr3s(annotation_utils::StructuralRegion::CDR3, alignment_info_, clone_set_);
        seqan::SeqFileOut out(output_config_.cdr3_compressed_fasta.c_str());
        std::vector<seqan::CharString> headers;
        std::vector<seqan::Dna5String> regions;
        for(auto it = compressed_cdr3s.cbegin(); it != compressed_cdr3s.cend(); it++) {
            headers.push_back(GetCompressedRegionFname(annotation_utils::StructuralRegion::CDR3,
                                                       it->first, it->second));
            regions.push_back(it->first.cdr_seq);
        }
        seqan::writeRecords(out, headers, regions);
        INFO("Compressed " << annotation_utils::StructuralRegion::CDR3 << " were written to " <<
                     output_config_.cdr3_compressed_fasta);
    }
}