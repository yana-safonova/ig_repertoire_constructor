#include <string>
#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <segfault_handler.hpp>
#include <logger/log_writers.hpp>
#include <seqan/seq_io.h>
#include "../umi_experiments/utils.hpp"
#include "../umi_experiments/umi_utils.hpp"
#include "../umi_experiments/utils/io.hpp"
#include "../fast_ig_tools/ig_matcher.hpp"
#include "ion_utils.hpp"

namespace {
    struct Params {
        std::string vjf_output_path;
        std::string alignment_info_path;
        std::string original_reads_path;
        std::string v_gene_path;
        std::string j_gene_path;
        std::string output_path;
    };

    bool read_args(int argc, const char *const *argv, Params& params) {
        namespace po = boost::program_options;
        po::options_description cmdl_options;
        cmdl_options.add_options()
                ("help,h", "print help message")
                ("cleaned,c", po::value<std::string>(&params.vjf_output_path)->required(), "file with reads cleaned by VJFinder")
                ("align,a", po::value<std::string>(&params.alignment_info_path)->required(), "file with alignment info reported by VJFinder")
                ("reads,r", po::value<std::string>(&params.original_reads_path)->required(), "file with original uncompressed reads")
                ("v", po::value<std::string>(&params.v_gene_path)->required(), "file with V genes")
                ("j", po::value<std::string>(&params.j_gene_path)->required(), "file with J genes")
                ("output,o", po::value<std::string>(&params.output_path)->required(), "output file with uncompressed reads");
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(cmdl_options).run(), vm);
        if (vm.count("help") || argc == 1) {
            std::cout << cmdl_options << std::endl;
            return false;
        }
        po::notify(vm);
        return true;
    }

    struct Input {
        std::vector<seqan::CharString> cleaned_ids;
        std::vector<seqan::CharString> cleaned_ids_with_umis;
        std::vector<seqan::Dna5String> cleaned_reads;
        std::unordered_map<seqan::CharString, seqan::Dna5String> original_reads;
        std::unordered_map<seqan::CharString, seqan::Dna5String> germline_v_genes;
        std::unordered_map<seqan::CharString, seqan::Dna5String> germline_j_genes;
        std::vector<seqan::CharString> read_v_gene;
        std::vector<seqan::CharString> read_j_gene;
    };

    void parse_alignment_info(const std::string& input_file_name,
                              std::vector<seqan::CharString>& v_gene, std::vector<seqan::CharString>& j_gene) {
        std::ifstream alignment(input_file_name);
        std::string line;
        std::getline(alignment, line);
        const auto columns = split(line, '\t');
        size_t v_gene_pos = (size_t) (std::find(columns.begin(), columns.end(), "V_hit") - columns.begin());
        size_t j_gene_pos = (size_t) (std::find(columns.begin(), columns.end(), "J_hit") - columns.begin());
        VERIFY(v_gene_pos < columns.size() && j_gene_pos < columns.size());
        while (std::getline(alignment, line)) {
            const auto values = split(line, '\t');
            VERIFY(values.size() == columns.size());
            v_gene.emplace_back(values[v_gene_pos]);
            j_gene.emplace_back(values[j_gene_pos]);
            VERIFY(length(*v_gene.crbegin()) > 0);
            VERIFY(length(*j_gene.crbegin()) > 0);
        }
    }

    void read_input(Params& params, Input& input) {
        read_seqan_records(params.vjf_output_path, input.cleaned_ids_with_umis, input.cleaned_reads);
        input.cleaned_reads.clear();
        read_seqan_records_cut_umi(params.vjf_output_path, input.cleaned_ids, input.cleaned_reads);
        parse_alignment_info(params.alignment_info_path, input.read_v_gene, input.read_j_gene);
        VERIFY(input.cleaned_reads.size() == input.read_v_gene.size());
        read_seqan_records_cut_umi(params.original_reads_path, input.original_reads);
        read_seqan_records_cut_umi(params.v_gene_path, input.germline_v_genes);
        read_seqan_records_cut_umi(params.j_gene_path, input.germline_j_genes);
    }

    std::unordered_map<size_t, size_t> hash_to_pos(const seqan::Dna5String& s) {
        const size_t K = 30;
        const auto comp_hashes = polyhashes(s, K);
        std::unordered_map<size_t, size_t> map;
        for (size_t i = 0; i < comp_hashes.size(); i ++) {
            VERIFY(map.find(comp_hashes[i]) == map.end());
            map[comp_hashes[i]] = i;
        }
        return map;
    };

    void skip_equal(const seqan::Dna5String& s, size_t& pos, boost::optional<std::stringstream&> out = {}) {
        const auto current = s[pos];
        while (pos < length(s) && s[pos] == current) {
            if (out) {
                (*out) << current;
            }
            pos ++;
        }
    }

    void cut_by_compressed_indices(const seqan::Dna5String& s, const size_t from, const size_t to, std::stringstream& out) {
        size_t pos = 0;
        for (size_t i = 0; i < from; i ++) {
            VERIFY_MSG(pos < length(s), "from " << from);
            skip_equal(s, pos);
        }
        for (size_t i = from; i < to; i ++) {
            VERIFY_MSG(pos < length(s), "to " << to);
            skip_equal(s, pos, { out });
        }
    }

    void get_common_kmers(const std::unordered_map<size_t, size_t>& comp_hash_to_pos,
            const std::unordered_map<size_t, size_t>& orig_comp_hash_to_pos, size_t& comp_pos, size_t& orig_pos) {
        for (const auto& entry : comp_hash_to_pos) {
            auto hash = entry.first;
            if (orig_comp_hash_to_pos.find(hash) != orig_comp_hash_to_pos.end()) {
                comp_pos = entry.second;
                orig_pos = orig_comp_hash_to_pos.at(hash);
                break;
            }
        }
    }
    
    seqan::Dna5String uncompress(const seqan::Dna5String& compressed, const seqan::Dna5String& orig,
                                 const seqan::Dna5String& v_gene, const seqan::Dna5String& j_gene) {
        // * Search for a pair of equal 30-mers in compressed original read and vjf cleaned read (also compressed) trying straight and reverse complement.
        // * Expand them by equality in both directions.
        // * Fill extents from germline.
        // * Uncompress.
        // All failures produce fatal error.

        // BTW, vjf output may not always be compressed and may contain repeating chars. But it seems to be OK.
        auto original = orig;
        VERIFY(length(original) > 0);
        const auto comp_hash_to_pos = hash_to_pos(compressed);
        auto original_compressed = compress_string(original);
        auto orig_comp_hash_to_pos = hash_to_pos(original_compressed);
        const size_t INF = std::numeric_limits<size_t>::max();
        size_t comp_pos = INF;
        size_t orig_pos = INF;
        get_common_kmers(comp_hash_to_pos, orig_comp_hash_to_pos, comp_pos, orig_pos);
        if (comp_pos == INF || orig_pos == INF) {
            seqan::reverseComplement(original_compressed);
            seqan::reverseComplement(original);
            orig_comp_hash_to_pos = hash_to_pos(original_compressed);
            get_common_kmers(comp_hash_to_pos, orig_comp_hash_to_pos, comp_pos, orig_pos);
        }
        VERIFY(comp_pos != INF && orig_pos != INF);
        size_t comp_left = comp_pos;
        size_t comp_right = comp_pos;
        size_t orig_left = orig_pos;
        size_t orig_right = orig_pos;
        while (comp_left > 0 && orig_left > 0 && compressed[comp_left - 1] == original_compressed[orig_left - 1]) {
            comp_left --;
            orig_left --;
        }
        while (comp_right + 1 < length(compressed) && orig_right + 1 < length(original_compressed) &&
                compressed[comp_right + 1] == original_compressed[orig_right + 1]) {
            comp_right ++;
            orig_right ++;
        }
        std::stringstream out;
        cut_by_compressed_indices(v_gene, 0, comp_left, out);
        cut_by_compressed_indices(original, orig_left, orig_right, out);
        cut_by_compressed_indices(j_gene, length(j_gene) + comp_right + 1 - length(compressed), length(j_gene), out);
        return out.str();
    }
}

int main(int argc, const char* const* argv) {
    segfault_handler sh;
    create_console_logger(logging::L_TRACE);

    Params params;
    if (!read_args(argc, argv, params)) {
        return 1;
    }

    Input input;
    read_input(params, input);

    std::vector<seqan::Dna5String> uncompressed_reads;
    uncompressed_reads.reserve(input.cleaned_ids.size());
    for (size_t i = 0; i < input.cleaned_ids.size(); i ++) {
        VERIFY(input.original_reads.find(input.cleaned_ids[i]) != input.original_reads.end());
        uncompressed_reads.push_back(uncompress(input.cleaned_reads[i], input.original_reads[input.cleaned_ids[i]],
                                                input.germline_v_genes[input.read_v_gene[i]],
                                                input.germline_j_genes[input.read_j_gene[i]]));
    }

    INFO("Writing uncompressed records to " << params.output_path);
    write_seqan_records(params.output_path, input.cleaned_ids_with_umis, uncompressed_reads);

    return 0;
}

