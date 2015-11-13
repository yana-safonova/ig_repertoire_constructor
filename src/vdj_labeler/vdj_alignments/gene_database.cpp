#include "gene_database.hpp"

std::string IgGeneTypeToString(IgGeneType gene_type) {
    if(gene_type == variable_gene)
        return "variable";
    if(gene_type == diversity_gene)
        return "diversity";
    return "join";
}

//-------------------------------------------------------------------------------------------

std::ostream& operator<< (std::ostream &out, const IgGene &obj) {
    out << "Name: " << obj.name() << ". Seq (len " << obj.length() << "): " << obj.seq();
    return out;
}

//-------------------------------------------------------------------------------------------

void IgGeneDatabase::AddGenesFromFile(std::string filename) {
    std::vector<CharString> read_headers;
    std::vector<Dna5String> read_seqs;
    seqan::SeqFileIn seqFileIn_reads(filename.c_str());
    readRecords(read_headers, read_seqs, seqFileIn_reads);
    for(size_t i = 0; i < read_headers.size(); i++) {
        ig_genes_.push_back(IgGene(read_headers[i], read_seqs[i]));
    }
}

const IgGene& IgGeneDatabase::GetByIndex(size_t index) const {
    assert(index < size());
    return ig_genes_[index];
}

const IgGene& IgGeneDatabase::GetByName(std::string gene_name) const {
    return GetByName(CharString(gene_name.c_str()));
}

const IgGene& IgGeneDatabase::GetByName(CharString gene_name) const {
    for(auto it = cbegin(); it != cend(); it++)
        if(it->name() == gene_name)
            return *it;
    return IgGene();
}

std::ostream& operator<<(std::ostream &out, const IgGeneDatabase &ig_gene_db) {
    out << "Ig genes database. Gene type: " << IgGeneTypeToString(ig_gene_db.GeneType()) << ". # records: " <<
            ig_gene_db.size() << std::endl;
    for(auto it = ig_gene_db.cbegin(); it != ig_gene_db.cend(); it++)
        out << *it << std::endl;
    return out;
}

//-------------------------------------------------------------------------------------------

void HC_GenesDatabase::AddGenesFromFile(IgGeneType gene_type, std::string filename){
    if(gene_type == variable_gene)
        variable_genes_.AddGenesFromFile(filename);
    else
    if(gene_type == diversity_gene)
        diversity_genes_.AddGenesFromFile(filename);
    else
        join_genes_.AddGenesFromFile(filename);
}

size_t HC_GenesDatabase::GenesNumber(IgGeneType gene_type) const {
    if(gene_type == variable_gene)
        return variable_genes_.size();
    if(gene_type == diversity_gene)
        return diversity_genes_.size();
    return join_genes_.size();
}

IgGene HC_GenesDatabase::GetByIndex(IgGeneType gene_type, size_t index) const {
    if(gene_type == variable_gene)
        return variable_genes_.GetByIndex(index);
    if(gene_type == diversity_gene)
        return diversity_genes_.GetByIndex(index);
    return join_genes_.GetByIndex(index);
}

std::ostream& operator<<(std::ostream &out, const HC_GenesDatabase& obj) {
    out << obj.VariableGenes();
    out << "--------------------" << std::endl;
    out << obj.DiversityGenes();
    out << "--------------------" << std::endl;
    out << obj.JoinGenes();
    return out;
}

//-------------------------------------------------------------------------------------------

void LC_GenesDatabase::AddGenesFromFile(IgGeneType gene_type, std::string filename){
    if(gene_type == variable_gene)
        variable_genes_.AddGenesFromFile(filename);
    else
        join_genes_.AddGenesFromFile(filename);
}

size_t LC_GenesDatabase::GenesNumber(IgGeneType gene_type) const {
    if(gene_type == variable_gene)
        return variable_genes_.size();
    return join_genes_.size();
}

IgGene LC_GenesDatabase::GetByIndex(IgGeneType gene_type, size_t index) const {
    if(gene_type == variable_gene)
        return variable_genes_.GetByIndex(index);
    return join_genes_.GetByIndex(index);
}

std::ostream& operator<<(std::ostream &out, const LC_GenesDatabase& obj) {
    out << obj.VariableGenes();
    out << "--------------------" << std::endl;
    out << obj.JoinGenes();
    return out;
}