#pragma once

#include "include_me.hpp"
#include "utils/fasta_reader.hpp"
#include "utils/string_tools.hpp"


enum IgGeneType {variable_gene, diversity_gene, join_gene};

string IgGeneTypeToString(IgGeneType gene_type) {
    if(gene_type == variable_gene)
        return "variable";
    if(gene_type == diversity_gene)
        return "diversity";
    return "join";
}

// ----------------------------------------------------------------------------

class IgGene {
    //IgGeneType gene_type_;
    string gene_name_;
    string short_gene_name_;
    string gene_seq_;

public:
    IgGene() : gene_name_(), gene_seq_() { }

    IgGene(string gene_name, string gene_seq) :
            gene_name_(gene_name),
            short_gene_name_(""),
            gene_seq_(gene_seq) { }

    IgGene(string gene_name, string short_gene_name, string gene_seq) :
            gene_name_(gene_name),
            short_gene_name_(short_gene_name),
            gene_seq_(gene_seq) { }

    string GeneName() const { return gene_name_; }

    string ShortGeneName() const { return short_gene_name_; }

    string GeneSeq() const { return gene_seq_; }
};

ostream& operator<< (ostream &out, const IgGene &obj) {
    out << "Full name: " << obj.GeneName() << endl;
    out << "Short name: " << obj.ShortGeneName() << ". Seq: " << obj.GeneSeq();
    return out;
}

// ----------------------------------------------------------------------------
class ShortGeneNameExtractor {
public:
    virtual string ExtractShortName(string long_name) const = 0;
};

class TrivialShortGeneNameExtractor : public ShortGeneNameExtractor {
public:
    TrivialShortGeneNameExtractor() { }

    string ExtractShortName(string long_name) const {
        return long_name;
    }
};

class IMGTShortGeneNameExtractor : public ShortGeneNameExtractor {
public:
    IMGTShortGeneNameExtractor() { }

    string ExtractShortName(string long_name) const {
        auto splits = split(long_name, '|');
        assert(splits.size() > 1);
        return splits[1];
    }
};

typedef shared_ptr<ShortGeneNameExtractor> ShortGeneNameExtractorPtr;

// ----------------------------------------------------------------------------
class IgGeneDatabase {
    IgGeneType gene_type_;
    ShortGeneNameExtractorPtr gene_name_extractor_ptr_;
    vector<IgGene> ig_genes_;

public:
    IgGeneDatabase(IgGeneType gene_type, ShortGeneNameExtractorPtr gene_name_extractor_ptr) :
            gene_type_(gene_type),
            gene_name_extractor_ptr_(gene_name_extractor_ptr) { }

    void AddGenesFromFile(string filename) {
        SingleFastaReader fasta_reader(filename);
        auto reads = fasta_reader.Read();
        for(auto read = reads.begin(); read != reads.end(); read++)
            ig_genes_.push_back(IgGene(read->name, gene_name_extractor_ptr_->ExtractShortName(read->name), read->seq));
    }

    void Print(ostream &out) const {
        out << "Ig genes database. Gene type: " << IgGeneTypeToString(gene_type_) <<
                ". # records: " << ig_genes_.size() << endl;
        for(auto it = ig_genes_.begin(); it != ig_genes_.end(); it++)
            out << *it << endl;
    }

    size_t GenesNumber() const { return ig_genes_.size(); }

    IgGene GetByIndex(size_t index) const {
        assert(index < GenesNumber());
        return ig_genes_[index];
    }
};

// ----------------------------------------------------------------------------
class HC_GenesDatabase {
    IgGeneDatabase variable_genes_;
    IgGeneDatabase diversity_genes_;
    IgGeneDatabase join_genes_;

public:
    HC_GenesDatabase(ShortGeneNameExtractorPtr short_gene_name_extractor_ptr):
        variable_genes_(variable_gene, short_gene_name_extractor_ptr),
        diversity_genes_(diversity_gene, short_gene_name_extractor_ptr),
        join_genes_(join_gene, short_gene_name_extractor_ptr) { }

    void AddGenesFromFile(IgGeneType gene_type, string filename){
        if(gene_type == variable_gene)
            variable_genes_.AddGenesFromFile(filename);
        else
            if(gene_type == diversity_gene)
                diversity_genes_.AddGenesFromFile(filename);
            else
                join_genes_.AddGenesFromFile(filename);
    }

    size_t GenesNumber(IgGeneType gene_type) const {
        if(gene_type == variable_gene)
            return variable_genes_.GenesNumber();
        if(gene_type == diversity_gene)
            return diversity_genes_.GenesNumber();
        return join_genes_.GenesNumber();
    }

    IgGene GetByIndex(IgGeneType gene_type, size_t index) const {
        if(gene_type == variable_gene)
            return variable_genes_.GetByIndex(index);
        if(gene_type == diversity_gene)
            return diversity_genes_.GetByIndex(index);
        return join_genes_.GetByIndex(index);
    }

    void Print(ostream &out) const {
        variable_genes_.Print(out);
        out << "--------------------" << endl;
        diversity_genes_.Print(out);
        out << "--------------------" << endl;
        join_genes_.Print(out);
    }
};

typedef shared_ptr<HC_GenesDatabase> HC_GenesDatabase_Ptr;

// ----------------------------------------------------------------------------
class LC_GenesDatabase {
    IgGeneDatabase variable_genes_;
    IgGeneDatabase join_genes_;

public:
    LC_GenesDatabase(ShortGeneNameExtractorPtr short_gene_name_extractor_ptr):
            variable_genes_(variable_gene, short_gene_name_extractor_ptr),
            join_genes_(join_gene, short_gene_name_extractor_ptr) { }

    void AddGenesFromFile(IgGeneType gene_type, string filename){
        if(gene_type == variable_gene)
            variable_genes_.AddGenesFromFile(filename);
        else
            join_genes_.AddGenesFromFile(filename);
    }

    size_t GenesNumber(IgGeneType gene_type) const {
        if(gene_type == variable_gene)
            return variable_genes_.GenesNumber();
        return join_genes_.GenesNumber();
    }

    IgGene GetByIndex(IgGeneType gene_type, size_t index) const {
        if(gene_type == variable_gene)
            return variable_genes_.GetByIndex(index);
        return join_genes_.GetByIndex(index);
    }

    void Print(ostream &out) const {
        variable_genes_.Print(out);
        out << "--------------------" << endl;
        join_genes_.Print(out);
    }
};

typedef shared_ptr<LC_GenesDatabase> LC_GenesDatabase_Ptr;
