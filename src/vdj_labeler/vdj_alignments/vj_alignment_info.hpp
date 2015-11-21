#pragma once

#include "gene_database.hpp"
#include "alignment_structs.hpp"
#include "../fastq_read_archive.hpp"

struct VJFinderMagicConsts {
    static const int read_name_index = 0;
    static const int first_v_match_pos_on_read = 3;
    static const int first_v_match_pos_on_gene = 4;
    static const int last_v_match_pos_on_read = 5;
    static const int last_v_match_pos_on_gene = 6;
    static const int v_gene_name_index = 8;
    static const int first_j_match_pos_on_read = 11;
    static const int first_j_match_pos_on_gene = 12;
    static const int last_j_match_pos_on_read = 13;
    static const int last_j_match_pos_on_gene = 14;
    static const int j_gene_name_index = 16;
};

class VJAlignmentInfo {
    const IgGeneDatabase& v_gene_db_;
    const IgGeneDatabase& j_gene_db_;
    const FastqReadArchive& read_archive_;

    std::vector<IgGeneAlignmentPositions> v_segments_;
    std::vector<IgGeneAlignmentPositions> j_segments_;

    ReadPtr GetReadName(const std::vector<std::string> &splits);

    void AddVAlignment(const std::vector<std::string> &splits, ReadPtr read);

    void AddJAlignment(const std::vector<std::string> &splits, ReadPtr read);

    void ParseLine(std::string line);

public:
    VJAlignmentInfo(const IgGeneDatabase& v_gene_db,
                    const IgGeneDatabase& j_gene_db,
                    const FastqReadArchive& read_archive) :
        v_gene_db_(v_gene_db),
        j_gene_db_(j_gene_db),
        read_archive_(read_archive) { }

    void ExtractAlignment(std::string filename);

    size_t size() const {
        return v_segments_.size();
    }

    IgGeneAlignmentPositions GetVAlignmentByIndex(size_t index) const;

    IgGeneAlignmentPositions GetJAlignmentByIndex(size_t index) const;
};

std::ostream& operator<<(std::ostream& out, const VJAlignmentInfo& obj);