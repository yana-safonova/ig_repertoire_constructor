#include <verify.hpp>
#include <logger/logger.hpp>

#include "immune_gene_database.hpp"

#include <seqan/translation.h>

namespace germline_utils {
    void ImmuneGene::ComputeAASeq() {
        using namespace seqan;
        StringSet<String<AminoAcid>, Owner<ConcatDirect<> > > aa_seqs;
        translate(aa_seqs, gene_seq_, WITH_FRAME_SHIFTS);
        aa_seq_ = aa_seqs[orf_];
    }

    void ImmuneGene::SetORF(unsigned orf) {
        VERIFY_MSG(orf == 0 or orf == 1 or orf == 2, "ORF value " << orf << " is not valid");
        orf_ = orf;
        ComputeAASeq();
    }

    std::ostream &operator<<(std::ostream &out, const ImmuneGene &obj) {
        out << "Name: " << obj.name() << ", type: " << obj.GeneType() << std::endl;
        out << "Seq (len " << obj.length() << "): " << obj.seq();
        return out;
    }

    void ImmuneGeneDatabase::CheckConsistencyFatal(const ImmuneGene &immune_gene) {
        VERIFY_MSG(immune_gene.GeneType() == gene_type_, "Type of gene is not " << gene_type_);
    }

    size_t ImmuneGeneDatabase::AddGenesFromFile(std::string filename) {
        std::vector <seqan::CharString> read_headers;
        std::vector <seqan::Dna5String> read_seqs;
        seqan::SeqFileIn seqFileIn_reads(filename.c_str());
        seqan::readRecords(read_headers, read_seqs, seqFileIn_reads);
        for (size_t i = 0; i < read_headers.size(); i++) {
            ImmuneGene immune_gene(gene_type_, read_headers[i], read_seqs[i], i);
            AddImmuneGene(immune_gene);
        }
        INFO(read_headers.size() << " records were extracted from " << filename);
        return read_headers.size();
    }

    void ImmuneGeneDatabase::AddImmuneGene(ImmuneGene immune_gene) {
        CheckConsistencyFatal(immune_gene);
        immune_genes_.push_back(immune_gene);
        gene_name_map_[std::string(seqan::toCString(immune_gene.name()))] = immune_genes_.size() - 1;
    }

    const ImmuneGene &ImmuneGeneDatabase::operator[](size_t index) const {
        VERIFY_MSG(index < size(), "Index " << index << " exceeds number of records in immune gene DB");
        return immune_genes_[index];
    }

    ImmuneGene& ImmuneGeneDatabase::GetImmuneGeneByIndex(size_t index) {
        VERIFY_MSG(index < size(), "Index " << index << " exceeds number of records in immune gene DB");
        return immune_genes_.at(index);
    }

    const ImmuneGene &ImmuneGeneDatabase::GetByName(std::string gene_name) const {
        VERIFY_MSG(gene_name_map_.find(gene_name) != gene_name_map_.end(), "Immune gene DB does not contain gene " << gene_name);
        return immune_genes_[GetIndexByName(gene_name)];
    }

    size_t ImmuneGeneDatabase::GetIndexByName(std::string gene_name) const {
        VERIFY_MSG(gene_name_map_.find(gene_name) != gene_name_map_.end(),
                   "Immune gene DB does not contain gene name " << gene_name);
        return gene_name_map_.at(gene_name);
    }

    size_t ImmuneGeneDatabase::GetIndexByName(seqan::CharString gene_name) const {
        return GetIndexByName(std::string(seqan::toCString(gene_name)));
    }

    std::ostream &operator<<(std::ostream &out, const ImmuneGeneDatabase &immune_gene_db) {
        out << "Immune genes database. Gene type: " << immune_gene_db.GeneType() <<
        ". # records: " << immune_gene_db.size() << std::endl;
        for (auto it = immune_gene_db.cbegin(); it != immune_gene_db.cend(); it++)
            out << *it << std::endl;
        return out;
    }

}
