#pragma once

#include <seqan/seq_io.h>
#include <unordered_map>

#include "../chain_type.hpp"
#include "../germline_gene_type.hpp"
#include <annotation_utils/aa_annotation/aa_annotation.hpp>

namespace germline_utils {

    class ImmuneGene {
        ImmuneGeneType gene_type_;
        seqan::CharString gene_name_;
        seqan::Dna5String gene_seq_;
        size_t id_;
        unsigned orf_;
        seqan::String<seqan::AminoAcid> aa_seq_;

        void ComputeAASeq();

    public:
        ImmuneGene() :
                gene_type_(),
                gene_name_(),
                gene_seq_(),
                id_(size_t(-1)),
                orf_(0),
                aa_seq_() { }

        ImmuneGene(ImmuneGeneType gene_type,
                   seqan::CharString gene_name,
                   seqan::Dna5String gene_seq,
                   size_t id) :
                gene_type_(gene_type),
                gene_name_(gene_name),
                gene_seq_(gene_seq),
                id_(id),
                orf_(0),
                aa_seq_() {
            ComputeAASeq();
        }

        void SetORF(unsigned orf);

        unsigned ORF() const { return orf_; }

        const seqan::CharString &name() const { return gene_name_; }

        const seqan::Dna5String &seq() const { return gene_seq_; }

        const seqan::String<seqan::AminoAcid>& aa_seq() const { return aa_seq_; }

        size_t length() const { return static_cast<size_t>(seqan::length(gene_seq_)); }

        ImmuneGeneType GeneType() const { return gene_type_; }

        ChainType Chain() const { return gene_type_.Chain(); }

        LymphocyteType Lymphocyte() const { return gene_type_.Lymphocyte(); }

        SegmentType Segment() const { return gene_type_.Segment(); }

        size_t id() const { return id_; }
    };

    typedef std::shared_ptr <ImmuneGene> ImmuneGenePtr;

    std::ostream &operator<<(std::ostream &out, const ImmuneGene &obj);

// ----------------------------------------------------------------------------

// class stores genes of fixed type, e.g., IGHV or TRBD
    class ImmuneGeneDatabase {
        ImmuneGeneType gene_type_;
        std::vector <ImmuneGene> immune_genes_;
        // map from gene_name to its index in immune_genes_ vector
        std::unordered_map <std::string, size_t> gene_name_map_;

        void CheckConsistencyFatal(const ImmuneGene &immune_gene);

    public:
        ImmuneGeneDatabase() : gene_type_() { }

        ImmuneGeneDatabase(ImmuneGeneType gene_type) :
                gene_type_(gene_type) { }

        // method returns number of added records
        size_t AddGenesFromFile(std::string filename);

        void AddImmuneGene(ImmuneGene immune_gene);

        ImmuneGeneType GeneType() const { return gene_type_; }

        LymphocyteType Lymphocyte() const { return gene_type_.Lymphocyte(); }

        ChainType Chain() const { return gene_type_.Chain(); }

        SegmentType Segment() const { return gene_type_.Segment(); }


        typedef std::vector<ImmuneGene>::const_iterator ImmuneGeneConstIterator;

        ImmuneGeneConstIterator cbegin() const { return immune_genes_.cbegin(); }

        ImmuneGeneConstIterator cend() const { return immune_genes_.cend(); }

        typedef std::vector<ImmuneGene>::iterator ImmuneGeneIterator;

        ImmuneGeneIterator begin() { return immune_genes_.begin(); }

        ImmuneGeneIterator end() { return immune_genes_.end(); }


        size_t size() const { return immune_genes_.size(); }

        const ImmuneGene &operator[](size_t index) const;

        ImmuneGene& GetImmuneGeneByIndex(size_t index);


        const ImmuneGene &GetByName(std::string gene_name) const;

        size_t GetIndexByName(std::string gene_name) const;

        size_t GetIndexByName(seqan::CharString gene_name) const;
    };

    std::ostream &operator<<(std::ostream &out, const ImmuneGeneDatabase &immune_gene_db);

}