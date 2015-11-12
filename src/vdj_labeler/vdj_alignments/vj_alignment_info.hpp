#pragma once

#include "gene_database.hpp"
#include "alignment_structs.hpp"

class VJAlignmentInfo {
    const IgGeneDatabase& v_gene_db_;
    const IgGeneDatabase& j_gene_db_;

    std::vector<IgGeneAlignment> v_segments_;
    std::vector<IgGeneAlignment> j_segments_;
    std::vector<std::string> read_names_;

    void AddReadName(const std::vector<std::string> &splits);

    void AddVAlignment(const std::vector<std::string> &splits);

    void AddJAlignment(const std::vector<std::string> &splits);

    void ParseLine(std::string line);

public:
    VJAlignmentInfo(const IgGeneDatabase& v_gene_db,
                    const IgGeneDatabase& j_gene_db) :
        v_gene_db_(v_gene_db),
        j_gene_db_(j_gene_db) { }

    void ExtractAlignment(std::string filename);

    size_t size() {
        return read_names_.size();
    }

    IgGeneAlignment GetVAlignmentByIndex(size_t index) const;

    IgGeneAlignment GetJAlignmentByIndex(size_t index) const;

    std::string GetReadNameByIndex(size_t index) const;
};

std::ostream& operator<<(std::ostream& out, VJAlignmentInfo& obj);