#pragma once
#include <fstream>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <antevolo_config.hpp>
#include "evolutionary_graph_utils/evolutionary_edge.hpp"

namespace antevolo {
    class shm_exception : public std::exception {
        std::string what_;
    public:
        shm_exception(std::string s) :
                what_(s) {}
        const char* what() const noexcept {
            return what_.c_str();
        }
    };

    class ShmModel {
        size_t k_;
        boost::unordered_map<std::string, std::vector<double>> mutation_probs_;
        const AntEvoloConfig::InputParams& config_;


        size_t nucl_and_kmer_to_index(char nucl, const std::string& kmer) ;

        double TableTransitionProbWithSubst(
                const seqan::Dna5String& seq,
                size_t central_pos,
                char nucl) ;

        double KmerTransitionProb(const std::string& from_kmer, const std::string& to_kmer) ;

        std::string GetKmerByCentralIndex(const seqan::Dna5String &seq, size_t central_pos) ;

        double SizeDependentTransitionProb (const seqan::Dna5String& src_seq,
                                            const seqan::Dna5String& dst_seq,
                                            size_t src_start_pos,
                                            size_t dst_start_pos,
                                            std::vector<size_t> diff_positions) ;

    public:
        ShmModel(size_t k, const AntEvoloConfig::InputParams& config)
                : k_(k),
                  config_(config) {
            std::ifstream inp(config_.model_input);
            std::string line;
            getline(inp, line);
            while(getline(inp, line)) {
                boost::trim(line);
                boost::tokenizer<boost::escaped_list_separator<char>> t(line);
                auto it = t.begin();
                std::string kmer = *it;
                mutation_probs_[kmer] = std::vector<double>(3);
                for (size_t i = 0; i < 3; ++i) {
                    it++;
                    mutation_probs_[kmer][i] = boost::lexical_cast<double>(*it);
                }
            }
        }
        
        double CDR3TransitionProb(const EvolutionaryEdge& edge) ;
    };
}