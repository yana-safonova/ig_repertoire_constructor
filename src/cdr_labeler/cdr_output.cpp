#include "cdr_output.hpp"

#include "compressed_cdr_set.hpp"

#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <annotation_utils/shm_filter.hpp>

namespace cdr_labeler {
    std::ostream& CDRLabelingWriter::OutputCloneRegion(std::ostream& out,
                                                       const annotation_utils::AnnotatedClone &clone,
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
        out << "Read_name\tChain_type\tV_hit\tJ_hit\tAA_seq\tHas_stop_codon\tIn-frame\tProductive\t"
                       "CDR1_nucls\tCDR1_start\tCDR1_end\tCDR2_nucls\tCDR2_start\tCDR2_end\t"
                       "CDR3_nucls\tCDR3_start\tCDR3_end" << std::endl;
        for(auto it = clone_set_.cbegin(); it != clone_set_.cend(); it++) {
            annotation_utils::AnnotatedClone cdr_clone = *it;
            auto vj_hit = alignment_info_.GetVJHitsByRead(cdr_clone.Read());
            out << cdr_clone.Read().name << "\t" << cdr_clone.ChainType() << "\t" <<
            it->VAlignment().subject().name() << "\t" <<
            it->JAlignment().subject().name() << "\t" << cdr_clone.AA() << "\t" <<
            cdr_clone.HasStopCodon() << "\t" << cdr_clone.InFrame() << "\t" << cdr_clone.Productive() << "\t";
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
        CompressedCDRSet compressed_cdr3s(annotation_utils::StructuralRegion::CDR3, clone_set_);
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

    void CDRLabelingWriter::OutputVGeneAlignment() const {
        std::ofstream out(output_config_.v_alignment_fasta);
        size_t index = 1;
        for(auto it = clone_set_.cbegin(); it != clone_set_.cend(); it++) {
            // todo: remove duplication with VJAlignmentInfoOutput
            auto subject_row = seqan::row(it->VAlignment().Alignment(), 0);
            auto query_row = seqan::row(it->VAlignment().Alignment(), 1);
            out << ">INDEX:" << index << "|READ:" << it->Read().name << "|START_POS:" <<
                    it->VAlignment().StartSubjectPosition() << "|END_POS:" <<
                    it->VAlignment().EndSubjectPosition() << std::endl;
            out << query_row << std::endl;
            out << ">INDEX:" << index << "|GENE:" << it->VAlignment().subject().name() <<
                    "|START_POS:" << it->VAlignment().StartQueryPosition() << "|END_POS:" <<
                    it->VAlignment().EndQueryPosition() << "|CHAIN_TYPE:" <<
                    it->VAlignment().subject().Chain() << std::endl;
            out << subject_row << std::endl;
            index++;
        }
        out.close();
        INFO("V alignments were written to " << output_config_.v_alignment_fasta);
    }

    void CDRLabelingWriter::OutputSHMsForRead(std::ostream& out, const annotation_utils::GeneSegmentSHMs &shms) const {
        size_t max_skipped_start = output_config_.shm_output_details.v_start_max_skipped;
        size_t max_skipped_end = output_config_.shm_output_details.v_end_max_skipped;
        if(shms.SegmentType() == germline_utils::SegmentType::JoinSegment) {
            max_skipped_start = output_config_.shm_output_details.j_start_max_skipped;
            max_skipped_end = output_config_.shm_output_details.j_end_max_skipped;
        }
        annotation_utils::PositionalSHMFilter shm_filter(shms, max_skipped_start, max_skipped_end);
        if(shm_filter.NumGoodSHMs() == 0)
            return;
        out << "Read_name:" << shms.Read().name << "\tRead_length:" << shms.Read().length() <<
                "\tGene_name:" << shms.ImmuneGene().name() << "\tGene_length:" << shms.ImmuneGene().length() <<
                "\tSegment:" << shms.SegmentType() << "\tChain_type:" << shms.ImmuneGene().Chain() << std::endl;
        for(auto it = shms.cbegin(); it != shms.cend(); it++)
            if(!shm_filter.FilterSHM(*it))
                out << it->shm_type << "\t" << it->read_nucl_pos << "\t" <<
                        it->gene_nucl_pos << "\t" << it->read_nucl << "\t" << it->gene_nucl << "\t" << it->read_aa <<
                        "\t" << it->gene_aa << "\t" << it->IsSynonymous() << "\t" << it->ToStopCodon() << std::endl;
    }

    void CDRLabelingWriter::OutputSHMs() const {
        std::ofstream out(output_config_.shm_output_details.shm_details);
        out << "SHM_type\tRead_pos\tGene_pos\tRead_nucl\tGene_nucl\tRead_aa\tGene_aa\tIs_synonymous\tTo_stop_codon\n";
        for(auto it = clone_set_.cbegin(); it != clone_set_.cend(); it++) {
            OutputSHMsForRead(out, it->VSHMs());
            OutputSHMsForRead(out, it->JSHMs());
        }
        out.close();
        INFO("SHM getails were written to " << output_config_.shm_output_details.shm_details);
    }
}