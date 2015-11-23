#include "standard.hpp"
#include "logger/log_writers.hpp"

#include "segfault_handler.hpp"
#include "stacktrace.hpp"
#include "memory_limit.hpp"
#include "copy_file.hpp"
#include "perfcounter.hpp"
#include "runtime_k.hpp"
#include "segfault_handler.hpp"

#include "fastq_read_archive.hpp"
#include "vdj_alignments/gene_database.hpp"
#include "vdj_alignments/vj_alignment_info.hpp"

#include "vdj_alignments/aligners/right_v_tail_aligner.hpp"
#include "vdj_alignments/aligners/left_j_tail_aligner.hpp"

#include "recombination_calculator/hc_model_based_recombination_calculator.hpp"

void create_console_logger() {
    using namespace logging;
    string log_props_file = "";
    logger *lg = create_logger(path::FileExists(log_props_file) ? log_props_file : "");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

int main(int, char**) {
    perf_counter pc;
    segfault_handler sh;
    create_console_logger();

    std::string fastq_reads_fname = "src/vdj_labeler/test/vdj_labeling.fastq";
    std::string vj_alignment_fname = "src/vdj_labeler/test/vdj_labeling.csv";
    std::string v_germline_genes_fname = "src/fast_ig_tools/germline/human/IGHV.fa";
    std::string d_germline_genes_fname = "src/fast_ig_tools/germline/human/IGHD.fa";
    std::string j_germline_genes_fname = "src/fast_ig_tools/germline/human/IGHJ.fa";

    INFO("VDJ labeler starts");

    FastqReadArchive reads_archive(fastq_reads_fname);
    INFO(reads_archive.size() << " reads were extracted from " << fastq_reads_fname);

    HC_GenesDatabase hc_db;
    hc_db.AddGenesFromFile(IgGeneType::variable_gene, v_germline_genes_fname);
    hc_db.AddGenesFromFile(IgGeneType::diversity_gene, d_germline_genes_fname);
    hc_db.AddGenesFromFile(IgGeneType::join_gene, j_germline_genes_fname);

    INFO(hc_db.VariableGenes().size() << " variable genes were extracted from " << v_germline_genes_fname);
    INFO(hc_db.DiversityGenes().size() << " diversity genes were extracted from " << d_germline_genes_fname);
    INFO(hc_db.JoinGenes().size() << " join genes were extracted from " << j_germline_genes_fname);

    VJAlignmentInfo vj_alignment_info(hc_db.VariableGenes(), hc_db.JoinGenes(), reads_archive);
    vj_alignment_info.ExtractAlignment(vj_alignment_fname);
    INFO(vj_alignment_info.size() << " alignment lines were extracted from " << vj_alignment_fname);
    INFO(vj_alignment_info);

    INFO("Alignment of right tails of V starts");
    RightVTailAligner raligner;
    for(size_t i = 0; i < vj_alignment_info.size(); i++) {
        auto v_alignment = raligner.ComputeAlignment(vj_alignment_info.GetVAlignmentByIndex(i));
        std::cout << *v_alignment << std::endl;
        std::cout << "---------" << std::endl;
    }
    INFO("Alignment of right tails of V ends");

    INFO("Alignment of left tails of J starts");
    LeftJTailAligner laligner;
    for(size_t i = 0; i < vj_alignment_info.size(); i++) {
        auto j_alignment = laligner.ComputeAlignment(vj_alignment_info.GetJAlignmentByIndex(i));
        std::cout << *j_alignment << std::endl;
        std::cout << "---------" << std::endl;
    }
    INFO("Alignment of left tails of J ends");

    INFO("Alignment of D segment starts");
    INFO("Alignment of D segment ends");

    INFO("VDJ labeler ends");
    unsigned ms = (unsigned)pc.time_ms();
    unsigned secs = (ms / 1000) % 60;
    unsigned mins = (ms / 1000 / 60) % 60;
    unsigned hours = (ms / 1000 / 60 / 60);
    INFO("Running time: " << hours << " hours " << mins << " minutes " << secs << " seconds");

    return 0;
}