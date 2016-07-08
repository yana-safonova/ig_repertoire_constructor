#pragma once

#include <seqan/seq_io.h>

#include "../chain_type.hpp"
#include "../germline_gene_type.hpp"

namespace germline_utils {

    class ImmuneGene {
        ImmuneGeneType gene_type_;
        seqan::CharString gene_name_;
        seqan::Dna5String gene_seq_;
        size_t id_;

    public:
        ImmuneGene() :
                gene_type_(),
                gene_name_(),
                gene_seq_(),
                id_(size_t(-1)) { }

        ImmuneGene(ImmuneGeneType gene_type,
                   seqan::CharString gene_name,
                   seqan::Dna5String gene_seq,
                   size_t id) :
                gene_type_(gene_type),
                gene_name_(gene_name),
                gene_seq_(gene_seq),
                id_(id) { }

        const seqan::CharString &name() const { return gene_name_; }

        const seqan::Dna5String &seq() const { return gene_seq_; }

        size_t length() const { return static_cast<size_t>(seqan::length(gene_seq_)); }

        ImmuneGeneType GeneType() const { return gene_type_; }

        size_t id() const { return id_; }
    };

    typedef std::shared_ptr <ImmuneGene> ImmuneGenePtr;

    std::ostream &operator<<(std::ostream &out, const ImmuneGene &obj);

// ----------------------------------------------------------------------------

// class stores genes of fixed type, e.g., IGHV or TRBD
    class ImmuneGeneDatabase {
        ImmuneGeneType gene_type_;
        std::vector <ImmuneGene> immune_genes_;
        // map from gene_name to its id
        std::unordered_map <std::string, size_t> gene_name_map_;

    public:
        ImmuneGeneDatabase(ImmuneGeneType gene_type) :
                gene_type_(gene_type) { }

        void AddGenesFromFile(std::string filename);


        ImmuneGeneType GeneType() const { return gene_type_; }

        LymphocyteType Lymphocyte() const { return gene_type_.Lymphocyte(); }

        ChainType Chain() const { return gene_type_.Chain(); }

        SegmentType Segment() const { return gene_type_.Segment(); }


        typedef std::vector<ImmuneGene>::const_iterator citerator;

        citerator cbegin() const { return ig_genes_.cbegin(); }

        citerator cend() const { return ig_genes_.cend(); }

        size_t size() const { return immune_genes_.size(); }

        const ImmuneGene &operator[](size_t index) const;


        const ImmuneGene &GetByName(std::string gene_name) const;

        size_t GetIndexByName(std::string gene_name) const;

        size_t GetIndexByName(seqan::CharString gene_name) const;
    };

    std::ostream &operator<<(std::ostream &out, const ImmuneGeneDatabase &immune_gene_db);

}