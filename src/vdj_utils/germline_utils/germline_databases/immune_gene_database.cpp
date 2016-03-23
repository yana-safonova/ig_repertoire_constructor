#include "immune_gene_database.hpp"

#include <verify.hpp>
#include "seqan/sequence.h"

namespace germline_utils {

    std::ostream &operator<<(std::ostream &out, const ImmuneGene &obj) {
        out << "Name: " << obj.name() << ", type: " << obj.GeneType() << std::endl;
        out << "Seq (len " << obj.length() << "): " << obj.seq();
        return out;
    }

    void ImmuneGeneDatabase::AddGenesFromFile(std::string filename) {
        std::vector <seqan::CharString> read_headers;
        std::vector <seqan::Dna5String> read_seqs;
        seqan::SeqFileIn seqFileIn_reads(filename.c_str());
        seqan::readRecords(read_headers, read_seqs, seqFileIn_reads);
        for (size_t i = 0; i < read_headers.size(); i++) {
            ImmuneGene immune_gene(gene_type_, read_headers[i], read_seqs[i], i);
            immune_genes_.push_back(immune_gene);
            gene_name_map_[std::string(seqan::toCString(read_headers[i]))] = i;
        }
    }

    const ImmuneGene &ImmuneGeneDatabase::operator[](size_t index) const {
        VERIFY_MSG(index < size(), "Index " << index << " exceeds number of records in immune gene DB");
        return immune_genes_[index];
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
