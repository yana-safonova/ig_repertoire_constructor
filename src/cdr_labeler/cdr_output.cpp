#include "cdr_output.hpp"

#include "compressed_cdr_set.hpp"
#include "CloneInfo.hpp"
#include <seqan/seq_io.h>
#include <seqan/stream.h>

namespace cdr_labeler {
    void CDRLabelingWriter::OutputCDRDetails() const {
        std::ofstream out(output_config_.feature_report_params.cdr_details);
        const auto columns = ReportColumns::ColumnSet<DivanReportEvalContext>::ChooseColumns(
                output_config_.feature_report_params.preset, output_config_.feature_report_params.columns
        );
        columns.PrintCsvHeader(out);
        size_t total_clone_sizes = 0;
        {
            for (const auto& cdr_clone : clone_set_) {
                const auto clone_info = CloneInfo::TryParse(cdr_clone.Read().name);
                if (clone_info) total_clone_sizes += clone_info->size;
            }
        }
        for (const auto& cdr_clone : clone_set_) {
            columns.Print(out, DivanReportEvalContext{cdr_clone, cdr_clone.VAlignment(), cdr_clone.JAlignment(), total_clone_sizes});
        }
        out.close();
        INFO("CDR details were written to " << output_config_.feature_report_params.cdr_details);
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
        //size_t max_skipped_start = output_config_.shm_output_details.v_start_max_skipped;
        //size_t max_skipped_end = output_config_.shm_output_details.v_end_max_skipped;
        //if(shms.SegmentType() == germline_utils::SegmentType::JoinSegment) {
        //    max_skipped_start = output_config_.shm_output_details.j_start_max_skipped;
        //    max_skipped_end = output_config_.shm_output_details.j_end_max_skipped;
        //}
        out << "Read_name:" << shms.Read().name << "\tRead_length:" << shms.Read().length() <<
                "\tGene_name:" << shms.ImmuneGene().name() << "\tGene_length:" << shms.ImmuneGene().length() <<
                "\tSegment:" << shms.SegmentType() << "\tChain_type:" << shms.ImmuneGene().Chain() << std::endl;
        for(auto it = shms.cbegin(); it != shms.cend(); it++) {
            out << it->shm_type << "\t" << it->read_nucl_pos << "\t" <<
            it->gene_nucl_pos << "\t" << it->read_nucl << "\t" << it->gene_nucl << "\t" << it->read_aa <<
            "\t" << it->gene_aa << "\t" << it->IsSynonymous() << "\t" << it->ToStopCodon() << std::endl;
        }
    }

    void CDRLabelingWriter::OutputSHMs() const {
        std::ofstream out(output_config_.shm_details);
        out << "SHM_type\tRead_pos\tGene_pos\tRead_nucl\tGene_nucl\tRead_aa\tGene_aa\tIs_synonymous\tTo_stop_codon\n";
        for(auto it = clone_set_.cbegin(); it != clone_set_.cend(); it++) {
            OutputSHMsForRead(out, it->VSHMs());
            OutputSHMsForRead(out, it->JSHMs());
        }
        out.close();
        INFO("SHM getails were written to " << output_config_.shm_details);
    }

    void CDRLabelingWriter::OutputCleanedReads() const {
        std::ofstream out(output_config_.cleaned_reads);
        for(auto it = clone_set_.cbegin(); it != clone_set_.cend(); it++) {
            out << ">" << it->Read().name << std::endl;
            out << it->Read().seq << std::endl;
        }
        out.close();
        INFO("Cleaned reads were written to " << output_config_.cleaned_reads);
    }
}

namespace ReportColumns {
    using DivanReportColumn = ReportColumns::Column<cdr_labeler::DivanReportEvalContext>;
    using DivanReportColumnSet = ReportColumns::ColumnSet<cdr_labeler::DivanReportEvalContext>;

    namespace DiversityAnalyzer {

        void print_region_nucleotides(std::ostream& out, const annotation_utils::AnnotatedClone& cdr_clone,
                const annotation_utils::StructuralRegion region) {
            if (cdr_clone.RegionIsEmpty(region)) {
                out << "-";
            } else {
                out << cdr_clone.GetRegionString(region);
            }
        }

        void print_region_aa(std::ostream& out, const annotation_utils::AnnotatedClone& cdr_clone,
                const annotation_utils::StructuralRegion region) {
            if (cdr_clone.RegionIsEmpty(region)) {
                out << "-";
            } else {
                seqan::String<seqan::AminoAcid> s;
                translate(s, cdr_clone.GetRegionString(region));
                out << s;
            }
        }

        void append_mixcr_alignment(std::ostream& out, const alignment_utils::ImmuneGeneReadAlignment& v_alignment,
                const annotation_utils::GeneSegmentSHMs& shms) {
            out << v_alignment.StartSubjectPosition() << '|'
                << v_alignment.EndSubjectPosition() + 1 << '|'
                << v_alignment.subject().length() << '|'
                << v_alignment.StartQueryPosition() << '|'
                << v_alignment.EndQueryPosition() + 1 << '|';
            shms.AppendInMixcrFormat(out);
            out << '|' << v_alignment.Score();
        }

        static const DivanReportColumn CLONE_NAME =
                {"Clone_name", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    out << context.cdr_clone.Read().name; }};
        static const DivanReportColumn CLONE_ID =
                {"Clone_ID", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    const auto clone_info = CloneInfo::TryParse(context.cdr_clone.Read().name);
                    VERIFY_MSG(clone_info, "Unknown clone name format. "
                            "Are you trying to use Clone_ID column for bare IgDiversityAnalyzer without running IgReC?");
                    out << clone_info->id;
                }};
        static const DivanReportColumn CLONE_COUNT =
                {"Clone_count", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    const auto clone_info = CloneInfo::TryParse(context.cdr_clone.Read().name);
                    VERIFY_MSG(clone_info, "Unknown clone name format. "
                            "Are you trying to use Clone_count column for bare IgDiversityAnalyzer without running IgReC?");
                    out << clone_info->size;
                }};
        static const DivanReportColumn CLONE_FRACTION =
                {"Clone_fraction", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    const auto clone_info = CloneInfo::TryParse(context.cdr_clone.Read().name);
                    VERIFY_MSG(clone_info, "Unknown clone name format. "
                            "Are you trying to use Clone_fraction column for bare IgDiversityAnalyzer without running IgReC?");
                    out << static_cast<double>(clone_info->size) / static_cast<double>(context.total_clone_sizes);
                }};
        static const DivanReportColumn CLONE_SEQUENCE =
                {"Clone_sequence", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    out << context.cdr_clone.Read().seq; }};
        static const DivanReportColumn CHAIN_TYPE =
                {"Chain_type", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    out << context.cdr_clone.ChainType(); }};
        static const DivanReportColumn V_HIT =
                {"V_hit", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    out << context.v_alignment.subject().name(); }};
        static const DivanReportColumn V_ALIGNMENT =
                {"V_alignment", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    append_mixcr_alignment(out, context.v_alignment, context.cdr_clone.VSHMs()); }};
        static const DivanReportColumn D_HIT =
                {"D_hit", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext&) {
                    out << "Unsupported"; }};
        static const DivanReportColumn D_ALIGNMENT =
                {"D_alignment", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext&) {
                    out << "Unsupported"; }};
        static const DivanReportColumn J_HIT =
                {"J_hit", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    out << context.j_alignment.subject().name(); }};
        static const DivanReportColumn J_ALIGNMENT =
                {"J_alignment", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    append_mixcr_alignment(out, context.j_alignment, context.cdr_clone.JSHMs()); }};
        static const DivanReportColumn C_HIT =
                {"C_hit", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext&) {
                    out << "Unsupported"; }};
        static const DivanReportColumn C_ALIGNMENT =
                {"C_alignment", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext&) {
                    out << "Unsupported"; }};
        static const DivanReportColumn AA_SEQ =
                {"AA_seq", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    out << context.cdr_clone.AA(); }};
        static const DivanReportColumn HAS_STOP_CODON =
                {"Has_stop_codon", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    out << context.cdr_clone.HasStopCodon(); }};
        static const DivanReportColumn IN_FRAME =
                {"In-frame", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    out << context.cdr_clone.InFrame(); }};
        static const DivanReportColumn PRODUCTIVE =
                {"Productive", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    out << context.cdr_clone.Productive(); }};
        static const DivanReportColumn FR1_NUCLS =
                {"FR1_nucls", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    print_region_nucleotides(out, context.cdr_clone, annotation_utils::StructuralRegion::FR1);
                }};
        static const DivanReportColumn FR1_AA =
                {"FR1_aa", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    print_region_aa(out, context.cdr_clone, annotation_utils::StructuralRegion::FR1);
                }};
        static const DivanReportColumn CDR1_NUCLS =
                {"CDR1_nucls", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    print_region_nucleotides(out, context.cdr_clone, annotation_utils::StructuralRegion::CDR1);
                }};
        static const DivanReportColumn CDR1_AA =
                {"CDR1_aa", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    print_region_aa(out, context.cdr_clone, annotation_utils::StructuralRegion::CDR1);
                }};
        static const DivanReportColumn CDR1_START =
                {"CDR1_start", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    const auto region = annotation_utils::StructuralRegion::CDR1;
                    if (context.cdr_clone.RegionIsEmpty(region)) {
                        out << "-";
                    } else {
                        out << context.cdr_clone.GetRangeByRegion(region).start_pos + 1;
                    }
                }};
        static const DivanReportColumn CDR1_END =
                {"CDR1_end", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    const auto region = annotation_utils::StructuralRegion::CDR1;
                    if (context.cdr_clone.RegionIsEmpty(region)) {
                        out << "-";
                    } else {
                        out << context.cdr_clone.GetRangeByRegion(region).end_pos + 1;
                    }
                }};
        static const DivanReportColumn FR2_NUCLS =
                {"FR2_nucls", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    print_region_nucleotides(out, context.cdr_clone, annotation_utils::StructuralRegion::FR2);
                }};
        static const DivanReportColumn FR2_AA =
                {"FR2_aa", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    print_region_aa(out, context.cdr_clone, annotation_utils::StructuralRegion::FR2);
                }};
        static const DivanReportColumn CDR2_NUCLS =
                {"CDR2_nucls", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    print_region_nucleotides(out, context.cdr_clone, annotation_utils::StructuralRegion::CDR2);
                }};
        static const DivanReportColumn CDR2_AA =
                {"CDR2_aa", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    print_region_aa(out, context.cdr_clone, annotation_utils::StructuralRegion::CDR2);
                }};
        static const DivanReportColumn CDR2_START =
                {"CDR2_start", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    const auto region = annotation_utils::StructuralRegion::CDR2;
                    if (context.cdr_clone.RegionIsEmpty(region)) {
                        out << "-";
                    } else {
                        out << context.cdr_clone.GetRangeByRegion(region).start_pos + 1;
                    }
                }};
        static const DivanReportColumn CDR2_END =
                {"CDR2_end", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    const auto region = annotation_utils::StructuralRegion::CDR2;
                    if (context.cdr_clone.RegionIsEmpty(region)) {
                        out << "-";
                    } else {
                        out << context.cdr_clone.GetRangeByRegion(region).end_pos + 1;
                    }
                }};
        static const DivanReportColumn FR3_NUCLS =
                {"FR3_nucls", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    print_region_nucleotides(out, context.cdr_clone, annotation_utils::StructuralRegion::FR3);
                }};
        static const DivanReportColumn FR3_AA =
                {"FR3_aa", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    print_region_aa(out, context.cdr_clone, annotation_utils::StructuralRegion::FR3);
                }};
        static const DivanReportColumn CDR3_NUCLS =
                {"CDR3_nucls", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    print_region_nucleotides(out, context.cdr_clone, annotation_utils::StructuralRegion::CDR3);
                }};
        static const DivanReportColumn CDR3_AA =
                {"CDR3_aa", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    print_region_aa(out, context.cdr_clone, annotation_utils::StructuralRegion::CDR3);
                }};
        static const DivanReportColumn CDR3_START =
                {"CDR3_start", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    const auto region = annotation_utils::StructuralRegion::CDR3;
                    if (context.cdr_clone.RegionIsEmpty(region)) {
                        out << "-";
                    } else {
                        out << context.cdr_clone.GetRangeByRegion(region).start_pos + 1;
                    }
                }};
        static const DivanReportColumn CDR3_END =
                {"CDR3_end", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    const auto region = annotation_utils::StructuralRegion::CDR3;
                    if (context.cdr_clone.RegionIsEmpty(region)) {
                        out << "-";
                    } else {
                        out << context.cdr_clone.GetRangeByRegion(region).end_pos + 1;
                    }
                }};
        static const DivanReportColumn FR4_NUCLS =
                {"FR4_nucls", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    print_region_nucleotides(out, context.cdr_clone, annotation_utils::StructuralRegion::FR4);
                }};
        static const DivanReportColumn FR4_AA =
                {"FR4_aa", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    print_region_aa(out, context.cdr_clone, annotation_utils::StructuralRegion::FR4);
                }};
        static const DivanReportColumn REF_POINTS =
                {"Ref_points", [](std::ostream& out, const cdr_labeler::DivanReportEvalContext& context) {
                    out << "::::";
                    if (! context.cdr_clone.RegionIsEmpty(annotation_utils::StructuralRegion::FR1)) {
                        out << context.cdr_clone.GetRangeByRegion(annotation_utils::StructuralRegion::FR1).start_pos;
                    }
                    out << ":";
                    if (! context.cdr_clone.RegionIsEmpty(annotation_utils::StructuralRegion::CDR1)) {
                        out << context.cdr_clone.GetRangeByRegion(annotation_utils::StructuralRegion::CDR1).start_pos;
                    }
                    out << ":";
                    if (! context.cdr_clone.RegionIsEmpty(annotation_utils::StructuralRegion::FR2)) {
                        out << context.cdr_clone.GetRangeByRegion(annotation_utils::StructuralRegion::FR2).start_pos;
                    }
                    out << ":";
                    if (! context.cdr_clone.RegionIsEmpty(annotation_utils::StructuralRegion::CDR2)) {
                        out << context.cdr_clone.GetRangeByRegion(annotation_utils::StructuralRegion::CDR2).start_pos;
                    }
                    out << ":";
                    if (! context.cdr_clone.RegionIsEmpty(annotation_utils::StructuralRegion::FR3)) {
                        out << context.cdr_clone.GetRangeByRegion(annotation_utils::StructuralRegion::FR3).start_pos;
                    }
                    out << ":";
                    if (! context.cdr_clone.RegionIsEmpty(annotation_utils::StructuralRegion::CDR3)) {
                        out << context.cdr_clone.GetRangeByRegion(annotation_utils::StructuralRegion::CDR3).start_pos;
                    }
                    // V deletion, V end, D begin, D 5' and 3' deletions, D end, J begin, J deletions skipped
                    out << ":::::::::";
                    if (! context.cdr_clone.RegionIsEmpty(annotation_utils::StructuralRegion::FR4)) {
                        out << context.cdr_clone.GetRangeByRegion(annotation_utils::StructuralRegion::FR4).start_pos;
                    }
                    out << ":";
                    if (! context.cdr_clone.RegionIsEmpty(annotation_utils::StructuralRegion::FR4)) {
                        out << context.cdr_clone.GetRangeByRegion(annotation_utils::StructuralRegion::FR4).end_pos + 1;
                    }
                    // C end skipped
                    out << ":";
                }};


        static const DivanReportColumnSet DIVAN_PRESET = {
                "divan",
                {
                        ReportColumns::DiversityAnalyzer::CLONE_NAME,
                        ReportColumns::DiversityAnalyzer::CHAIN_TYPE,
                        ReportColumns::DiversityAnalyzer::V_HIT,
                        ReportColumns::DiversityAnalyzer::J_HIT,
                        ReportColumns::DiversityAnalyzer::AA_SEQ,
                        ReportColumns::DiversityAnalyzer::HAS_STOP_CODON,
                        ReportColumns::DiversityAnalyzer::IN_FRAME,
                        ReportColumns::DiversityAnalyzer::PRODUCTIVE,
                        ReportColumns::DiversityAnalyzer::CDR1_NUCLS,
                        ReportColumns::DiversityAnalyzer::CDR1_START,
                        ReportColumns::DiversityAnalyzer::CDR1_END,
                        ReportColumns::DiversityAnalyzer::CDR2_NUCLS,
                        ReportColumns::DiversityAnalyzer::CDR2_START,
                        ReportColumns::DiversityAnalyzer::CDR2_END,
                        ReportColumns::DiversityAnalyzer::CDR3_NUCLS,
                        ReportColumns::DiversityAnalyzer::CDR3_START,
                        ReportColumns::DiversityAnalyzer::CDR3_END
                }
        };

        static const DivanReportColumnSet MIN_PRESET = {
                "min",
                {
                        ReportColumns::DiversityAnalyzer::CLONE_COUNT,
                        ReportColumns::DiversityAnalyzer::V_HIT,
                        ReportColumns::DiversityAnalyzer::D_HIT,
                        ReportColumns::DiversityAnalyzer::J_HIT,
                        ReportColumns::DiversityAnalyzer::C_HIT,
                        ReportColumns::DiversityAnalyzer::CDR3_NUCLS,
                }
        };

        static const DivanReportColumnSet MIXCR_FULL_PRESET = {
                "mixcr-full",
                {
                        ReportColumns::DiversityAnalyzer::CLONE_ID,
                        ReportColumns::DiversityAnalyzer::CLONE_COUNT,
                        ReportColumns::DiversityAnalyzer::CLONE_FRACTION,
                        ReportColumns::DiversityAnalyzer::CLONE_SEQUENCE,
                        ReportColumns::DiversityAnalyzer::V_HIT,
                        ReportColumns::DiversityAnalyzer::D_HIT,
                        ReportColumns::DiversityAnalyzer::J_HIT,
                        ReportColumns::DiversityAnalyzer::C_HIT,
                        ReportColumns::DiversityAnalyzer::V_ALIGNMENT,
                        ReportColumns::DiversityAnalyzer::D_ALIGNMENT,
                        ReportColumns::DiversityAnalyzer::J_ALIGNMENT,
                        ReportColumns::DiversityAnalyzer::C_ALIGNMENT,
                        ReportColumns::DiversityAnalyzer::FR1_NUCLS,
                        ReportColumns::DiversityAnalyzer::CDR1_NUCLS,
                        ReportColumns::DiversityAnalyzer::FR2_NUCLS,
                        ReportColumns::DiversityAnalyzer::CDR2_NUCLS,
                        ReportColumns::DiversityAnalyzer::FR3_NUCLS,
                        ReportColumns::DiversityAnalyzer::CDR3_NUCLS,
                        ReportColumns::DiversityAnalyzer::FR4_NUCLS,
                        ReportColumns::DiversityAnalyzer::FR1_AA,
                        ReportColumns::DiversityAnalyzer::CDR1_AA,
                        ReportColumns::DiversityAnalyzer::FR2_AA,
                        ReportColumns::DiversityAnalyzer::CDR2_AA,
                        ReportColumns::DiversityAnalyzer::FR3_AA,
                        ReportColumns::DiversityAnalyzer::CDR3_AA,
                        ReportColumns::DiversityAnalyzer::FR4_AA,
                        ReportColumns::DiversityAnalyzer::REF_POINTS,
                }
        };
    };

    template <>
    const std::vector<DivanReportColumn> DivanReportColumn::COLUMN_TYPES = {
            ReportColumns::DiversityAnalyzer::CLONE_NAME,
            ReportColumns::DiversityAnalyzer::CLONE_ID,
            ReportColumns::DiversityAnalyzer::CLONE_COUNT,
            ReportColumns::DiversityAnalyzer::CLONE_FRACTION,
            ReportColumns::DiversityAnalyzer::CLONE_SEQUENCE,
            ReportColumns::DiversityAnalyzer::CHAIN_TYPE,
            ReportColumns::DiversityAnalyzer::V_HIT,
            ReportColumns::DiversityAnalyzer::V_ALIGNMENT,
            ReportColumns::DiversityAnalyzer::D_HIT,
            ReportColumns::DiversityAnalyzer::D_ALIGNMENT,
            ReportColumns::DiversityAnalyzer::J_HIT,
            ReportColumns::DiversityAnalyzer::J_ALIGNMENT,
            ReportColumns::DiversityAnalyzer::C_HIT,
            ReportColumns::DiversityAnalyzer::C_ALIGNMENT,
            ReportColumns::DiversityAnalyzer::AA_SEQ,
            ReportColumns::DiversityAnalyzer::HAS_STOP_CODON,
            ReportColumns::DiversityAnalyzer::IN_FRAME,
            ReportColumns::DiversityAnalyzer::PRODUCTIVE,
            ReportColumns::DiversityAnalyzer::FR1_NUCLS,
            ReportColumns::DiversityAnalyzer::FR1_AA,
            ReportColumns::DiversityAnalyzer::CDR1_NUCLS,
            ReportColumns::DiversityAnalyzer::CDR1_AA,
            ReportColumns::DiversityAnalyzer::CDR1_START,
            ReportColumns::DiversityAnalyzer::CDR1_END,
            ReportColumns::DiversityAnalyzer::FR2_NUCLS,
            ReportColumns::DiversityAnalyzer::FR2_AA,
            ReportColumns::DiversityAnalyzer::CDR2_NUCLS,
            ReportColumns::DiversityAnalyzer::CDR2_AA,
            ReportColumns::DiversityAnalyzer::CDR2_START,
            ReportColumns::DiversityAnalyzer::CDR2_END,
            ReportColumns::DiversityAnalyzer::FR3_NUCLS,
            ReportColumns::DiversityAnalyzer::FR3_AA,
            ReportColumns::DiversityAnalyzer::CDR3_NUCLS,
            ReportColumns::DiversityAnalyzer::CDR3_AA,
            ReportColumns::DiversityAnalyzer::CDR3_START,
            ReportColumns::DiversityAnalyzer::CDR3_END,
            ReportColumns::DiversityAnalyzer::FR4_NUCLS,
            ReportColumns::DiversityAnalyzer::FR4_AA,
            ReportColumns::DiversityAnalyzer::REF_POINTS,
    };

    template <>
    const std::vector<DivanReportColumnSet> DivanReportColumnSet::PRESETS = {
            ReportColumns::DiversityAnalyzer::DIVAN_PRESET,
            ReportColumns::DiversityAnalyzer::MIN_PRESET,
            ReportColumns::DiversityAnalyzer::MIXCR_FULL_PRESET,
    };
}
