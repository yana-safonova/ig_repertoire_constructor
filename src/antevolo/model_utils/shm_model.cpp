#include "shm_model.hpp"


namespace antevolo {
    size_t ShmModel::nucl_and_kmer_to_index(char nucl, const std::string& kmer)  {
        switch(kmer[2]) {
            case 'A' :
                switch(nucl) {
                    case 'C' : return 0;
                    case 'G' : return 1;
                    case 'T' : return 2;
                    default: return 3;
                }
            case 'C' :
                switch(nucl) {
                    case 'A' : return 0;
                    case 'G' : return 1;
                    case 'T' : return 2;
                    default: return 3;
                }
            case 'G' :
                switch(nucl) {
                    case 'A' : return 0;
                    case 'C' : return 1;
                    case 'T' : return 2;
                    default: return 3;
                }
            case 'T' :
                switch(nucl) {
                    case 'A' : return 0;
                    case 'C' : return 1;
                    case 'G' : return 2;
                    default: return 3;
                }
            default: return 3;
        }
    }

    double ShmModel::TableTransitionProbWithSubst(
            const seqan::Dna5String &seq,
            size_t central_pos,
            char nucl)  {
        std::string kmer;
        for (size_t i = 0; i < k_; ++i) {
            kmer.push_back(char(seq[central_pos - (k_ - 1) / 2 + i]));
        }
        if (mutation_probs_.find(kmer) == mutation_probs_.end()) {
            INFO("Table: no such a kmer:'" << kmer.c_str() << "'" );
        }
        return mutation_probs_[kmer][nucl_and_kmer_to_index(nucl, kmer)];
    }

    double ShmModel::KmerTransitionProb(const std::string &from_kmer, const std::string &to_kmer)  {
        if (mutation_probs_.find(from_kmer) == mutation_probs_.end()) {
            INFO("kTable: no such a kmer:'" << from_kmer.c_str() << "'" );
        }
        return mutation_probs_[from_kmer][nucl_and_kmer_to_index(to_kmer[(k_ - 1) / 2], from_kmer)];
    }

    std::string ShmModel::GetKmerByCentralIndex(const seqan::Dna5String &seq, size_t central_pos)  {
        std::string kmer;
        for (size_t pos = central_pos - (k_ - 1) / 2; pos <= central_pos + (k_ - 1) / 2; ++pos) {
            kmer.push_back(seq[pos]);
        }
        return kmer;
    }

    double ShmModel::SizeDependentTransitionProb(const seqan::Dna5String &src_seq,
                                                 const seqan::Dna5String &dst_seq,
                                                 size_t src_start_pos,
                                                 size_t dst_start_pos,
                                                 std::vector<size_t> diff_positions)  {
        //std::cout << "Size" << std::endl;
        switch (diff_positions.size()) {
            case 0:
                return 1.0;
            case 1 : {
                //std::cout << "1" << std::endl;
                return TableTransitionProbWithSubst(src_seq,
                                                    src_start_pos + diff_positions[0],
                                                    dst_seq[dst_start_pos + diff_positions[0]]);
            }
            case 2 : {
                if (diff_positions.at(1) - diff_positions.at(0) > (k_ - 1) / 2) {
                    //std::cout << "2_1" << std::endl;
                    return TableTransitionProbWithSubst(src_seq,
                                                        src_start_pos + diff_positions[0],
                                                        dst_seq[dst_start_pos + diff_positions[0]]) *
                           TableTransitionProbWithSubst(src_seq,
                                                        src_start_pos + diff_positions[1],
                                                        dst_seq[dst_start_pos + diff_positions[1]]);
                }
                double best_way = 0.0;
                std::vector<std::vector<size_t>> mutation_orders = {{0, 1},
                                                                    {1, 0}};
                //std::cout << "2_2" << std::endl;
                for (auto const& order : mutation_orders) {
                    double way = 1.0;
                    std::string from_kmer_1 = GetKmerByCentralIndex(src_seq,
                                                                    src_start_pos + diff_positions[order[0]]);
                    std::string to_kmer_1(from_kmer_1);
                    to_kmer_1[(k_ - 1) / 2] = dst_seq[dst_start_pos + diff_positions[order[0]]];
                    way *= KmerTransitionProb(from_kmer_1, to_kmer_1);
                    std::string to_kmer_2 = GetKmerByCentralIndex(dst_seq, dst_start_pos + diff_positions[order[1]]);
                    std::string from_kmer_2(to_kmer_2);
                    from_kmer_2[(k_ - 1) / 2 + diff_positions[order[0]] - diff_positions[order[1]]] =
                            src_seq[src_start_pos + diff_positions[order[0]]];
                    way *= KmerTransitionProb(from_kmer_2, to_kmer_2);
                    best_way = std::max(best_way, way);
                }
                return best_way;
            }
            case 3: {
                //std::cout << diff_positions.size() << std::endl;
                if (diff_positions[1] - diff_positions[0] > (k_ - 1) / 2) {
                    return SizeDependentTransitionProb(src_seq,
                                                       dst_seq,
                                                       src_start_pos,
                                                       dst_start_pos,
                                                       {diff_positions[0]}) *
                           SizeDependentTransitionProb(src_seq,
                                                       dst_seq,
                                                       src_start_pos,
                                                       dst_start_pos,
                                                       {diff_positions[1], diff_positions[2]});
                }
                if (diff_positions[2] - diff_positions[1] > (k_ - 1) / 2) {
                    return SizeDependentTransitionProb(src_seq,
                                                       dst_seq,
                                                       src_start_pos,
                                                       dst_start_pos,
                                                       {diff_positions[0], diff_positions[1]}) *
                           SizeDependentTransitionProb(src_seq,
                                                       dst_seq,
                                                       src_start_pos,
                                                       dst_start_pos,
                                                       {diff_positions[2]});
                }
                //double best_way = 0;
                /*std::vector<std::vector<size_t>> mutation_orders =
                        {{0, 1, 2},
                         {0, 2, 1},
                         {1, 0, 2},
                         {1, 2, 0},
                         {2, 0, 1},
                         {2, 1, 0}};*/
                return SizeDependentTransitionProb(src_seq,
                                                   dst_seq,
                                                   src_start_pos,
                                                   dst_start_pos,
                                                   {diff_positions[0], diff_positions[1]}) *
                       SizeDependentTransitionProb(src_seq,
                                                   dst_seq,
                                                   src_start_pos,
                                                   dst_start_pos,
                                                   {diff_positions[2]});
            }
            default: return 0;
        }
    }

    double ShmModel::CDR3TransitionProb(const EvolutionaryEdge& edge)  {
        if (edge.src_clone.get() == nullptr || edge.dst_clone.get() == nullptr)
        {
            INFO("Invalid edge endpoints pointers");
            return 0.00001;
        }
        if (edge.src_clone->RegionIsEmpty(annotation_utils::StructuralRegion::CDR3) ||
                edge.dst_clone->RegionIsEmpty(annotation_utils::StructuralRegion::CDR3)) {
            INFO("CDR3 region is empty");
            return 0.00001;
        }
        size_t src_start_pos = edge.src_clone->GetRangeByRegion(
                annotation_utils::StructuralRegion::CDR3).start_pos;
        size_t dst_start_pos = edge.dst_clone->GetRangeByRegion(
                annotation_utils::StructuralRegion::CDR3).start_pos;
        size_t length = edge.src_clone->GetRangeByRegion(
                annotation_utils::StructuralRegion::CDR3).length();
        if (length != edge.dst_clone->GetRangeByRegion(
                annotation_utils::StructuralRegion::CDR3).length()) {
            std::stringstream ss;
            ss << "nonequal CDR3 lengths: " << edge.src_clone_num << " " << edge.dst_clone_num;
            throw shm_exception(ss.str());
        }
        if (src_start_pos < 2 ||
                dst_start_pos < 2 ||
                src_start_pos + length >= edge.src_clone->Read().length() ||
                dst_start_pos + length >= edge.dst_clone->Read().length()) {
            throw shm_exception("bad CDR3 ranges: no 2bp spaces around");
        }

        auto const& src_seq = edge.src_clone->Read().seq;
        auto const& dst_seq = edge.dst_clone->Read().seq;
        std::vector<size_t> diff_positions;

        if (src_seq[src_start_pos - 2] == 'N' || dst_seq[dst_start_pos - 2] == 'N' ||
            src_seq[src_start_pos - 1] == 'N' || dst_seq[dst_start_pos - 1] == 'N' ||
            src_seq[src_start_pos + length] == 'N' || dst_seq[dst_start_pos + length] == 'N' ||
            src_seq[src_start_pos + length + 1] == 'N' || dst_seq[dst_start_pos + length + 1] == 'N') {
            return 0.00001;
        }

        for (size_t pos = 0; pos < length; ++pos) {
            if (src_seq[src_start_pos + pos] == 'N' || dst_seq[dst_start_pos + pos] == 'N') {
                return 0.00001;
            }
            if (src_seq[src_start_pos + pos] != dst_seq[dst_start_pos + pos]) {
                diff_positions.push_back(pos);
            }
        }
        return SizeDependentTransitionProb(src_seq,
                                           dst_seq,
                                           src_start_pos,
                                           dst_start_pos,
                                           diff_positions);
    }
}