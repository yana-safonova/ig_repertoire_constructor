#pragma once

#include <cassert>
#include <vector>
#include <stdexcept>

using std::vector;
using std::cout;
using std::cin;
using std::cerr;
using std::endl;

#include <boost/program_options.hpp>
#include "fast_ig_tools.hpp"
using path::make_dirs;

#include <seqan/seq_io.h>
using seqan::Dna5String;
using seqan::CharString;
using seqan::length;
using seqan::SeqFileIn;
using seqan::SeqFileOut;

#include "ig_block_alignment.hpp"

namespace fast_ig_tools {
using std::string;

struct GermlineLociVJDB {
    std::vector<Dna5String> v_reads;
    std::vector<CharString> v_ids;
    std::vector<Dna5String> j_reads;
    std::vector<CharString> j_ids;
};

class GermlineDB {
public:
    std::vector<GermlineLociVJDB> locus_databases;
    std::vector<std::string> locus_names;
    const std::string db_directory;
    const std::string organism;
    bool pseudogenes;
    std::vector<Dna5String> all_v_reads;
    std::vector<Dna5String> all_j_reads;
    std::vector<CharString> all_v_ids;
    std::vector<CharString> all_j_ids;
    std::vector<size_t> locus_index;



    static std::vector<std::string> expand_loci(const std::vector<std::string> &loci) {
        std::vector<std::string> IG = { "IGH", "IGK", "IGL" };
        std::vector<std::string> TR = { "TRA", "TRB", "TRD", "TRG" };
        std::vector<std::string> ALL(IG);
        ALL.insert(ALL.end(), TR.cbegin(), TR.cend());

        std::vector<std::string> result;
        for (const std::string &locus : loci) {
            if (locus == "IG") {
                result.insert(result.end(), IG.cbegin(), IG.cend());
            } else if (locus == "TR") {
                result.insert(result.end(), TR.cbegin(), TR.cend());
            } else if (locus == "all") {
                result.insert(result.end(), ALL.cbegin(), ALL.cend());
            } else if (std::find(ALL.cbegin(), ALL.cend(), locus) != ALL.cend()) {
                result.push_back(locus);
            } else {
                throw std::invalid_argument("Invalid locus name: " + locus);
            }
        }

        // Remove duplicates
        remove_duplicates(result);

        return result;
    }

    std::string gene_file_name(const std::string &locus,
                               const std::string &gene,
                               bool pseudo = false) const {
        return db_directory + "/" + organism + "/" + locus + gene + (pseudo ? "-allP" : "") + ".fa";
    }

    GermlineDB(const std::string &db_directory,
               const std::string &organism,
               bool pseudogenes,
               const std::vector<std::string> &loci) : db_directory{db_directory},
                                                       organism{organism},
                                                       pseudogenes{pseudogenes},
                                                       locus_names{expand_loci(loci)} {
        for (const std::string &locus : locus_names) {
            GermlineLociVJDB db;

            {
                std::string v_file = gene_file_name(locus, "V", pseudogenes);
                SeqFileIn seqFileIn(v_file.c_str());
                readRecords(db.v_ids, db.v_reads, seqFileIn);
            }

            {
                std::string j_file = gene_file_name(locus, "J", pseudogenes);
                SeqFileIn seqFileIn(j_file.c_str());
                readRecords(db.j_ids, db.j_reads, seqFileIn);
            }

            locus_databases.push_back(db);
        }

        // Join V genes
        for (size_t i = 0; i < locus_databases.size(); ++i) {
            const auto &db = locus_databases[i];

            all_v_reads.insert(all_v_reads.end(),
                               db.v_reads.cbegin(), db.v_reads.cend());
            all_v_ids.insert(all_v_ids.end(),
                               db.v_ids.cbegin(), db.v_ids.cend());

            all_j_reads.insert(all_j_reads.end(),
                               db.j_reads.cbegin(), db.j_reads.cend());
            all_j_ids.insert(all_j_ids.end(),
                               db.j_ids.cbegin(), db.j_ids.cend());

            for (const auto &_ : db.v_reads) {
                locus_index.push_back(i);
            }
        }

    }
};



} // namespace fast_ig_tools

// vim: ts=4:sw=4
