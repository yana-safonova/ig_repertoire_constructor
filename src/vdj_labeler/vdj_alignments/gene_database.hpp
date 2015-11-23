#pragma once

#include <seqan/seq_io.h>

using seqan::Dna5String;
using seqan::CharString;
using seqan::length;
using seqan::SeqFileIn;
using seqan::SeqFileOut;

enum IgGeneType {variable_gene, diversity_gene, join_gene};

std::string IgGeneTypeToString(IgGeneType gene_type);

// ----------------------------------------------------------------------------

class IgGene {
    CharString gene_name_;
    Dna5String gene_seq_;

public:
    IgGene() : gene_name_(), gene_seq_() { }

    IgGene(CharString gene_name, Dna5String gene_seq) :
            gene_name_(gene_name),
            gene_seq_(gene_seq) { }

    CharString name() const { return gene_name_; }

    Dna5String seq() const { return gene_seq_; }

    size_t length() const { return static_cast<size_t>(seqan::length(gene_seq_)); }
};

typedef std::shared_ptr<IgGene> IgGenePtr;

std::ostream& operator<<(std::ostream &out, const IgGene &obj);

// ----------------------------------------------------------------------------

class IgGeneDatabase {
    IgGeneType gene_type_;
    std::vector<IgGenePtr> ig_genes_;
    std::map<std::string, IgGenePtr> gene_name_map_;

public:
    IgGeneDatabase(IgGeneType gene_type) :
            gene_type_(gene_type) { }

    void AddGenesFromFile(std::string filename);

    size_t size() const { return ig_genes_.size(); }

    IgGenePtr GetByIndex(size_t index) const;

    IgGeneType GeneType() const { return gene_type_; }

    typedef std::vector<IgGenePtr>::const_iterator citerator;

    citerator cbegin() const { return ig_genes_.cbegin(); }

    citerator cend() const { return ig_genes_.cend(); }

    IgGenePtr GetByName(std::string gene_name) const;
};

std::ostream& operator<<(std::ostream &out, const IgGeneDatabase &ig_gene_db);

// ----------------------------------------------------------------------------

class HC_GenesDatabase {
    IgGeneDatabase variable_genes_;
    IgGeneDatabase diversity_genes_;
    IgGeneDatabase join_genes_;

public:
    HC_GenesDatabase():
        variable_genes_(variable_gene),
        diversity_genes_(diversity_gene),
        join_genes_(join_gene) { }

    void AddGenesFromFile(IgGeneType gene_type, std::string filename);

    size_t GenesNumber(IgGeneType gene_type) const;

    IgGenePtr GetByIndex(IgGeneType gene_type, size_t index) const;

    const IgGeneDatabase& VariableGenes() const { return variable_genes_; }

    const IgGeneDatabase& DiversityGenes() const { return diversity_genes_; }

    const IgGeneDatabase& JoinGenes() const { return join_genes_; }
};

std::ostream& operator<<(std::ostream &out, const HC_GenesDatabase& obj);

typedef std::shared_ptr<HC_GenesDatabase> HC_GenesDatabase_Ptr;

// ----------------------------------------------------------------------------

class LC_GenesDatabase {
    IgGeneDatabase variable_genes_;
    IgGeneDatabase join_genes_;

public:
    LC_GenesDatabase():
            variable_genes_(variable_gene),
            join_genes_(join_gene) { }

    void AddGenesFromFile(IgGeneType gene_type, std::string filename);

    size_t GenesNumber(IgGeneType gene_type) const;

    IgGenePtr GetByIndex(IgGeneType gene_type, size_t index) const;

    const IgGeneDatabase& VariableGenes() const { return variable_genes_; }

    const IgGeneDatabase& JoinGenes() const { return join_genes_; }
};

std::ostream& operator<<(std::ostream &out, const LC_GenesDatabase& obj);

typedef std::shared_ptr<LC_GenesDatabase> LC_GenesDatabase_Ptr;