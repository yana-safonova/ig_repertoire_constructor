#include "parent_read_reconstructor.hpp"

namespace antevolo {

    static core::Read ParentReadReconstructor::ReconstructParentRead(
            const std::shared_ptr<annotation_utils::AnnotatedClone> &clone1,
            const std::shared_ptr<annotation_utils::AnnotatedClone> &clone2,
            size_t id) const {
        // read
        std::stringstream s;
        s << std::string("fake_read_") << id << ("___size___1");
        std::string read_name = s.str();
        std::string res_string;
        seqan::DnaString res_dna5string;

        auto read1 = clone1->Read().seq;
        auto read2 = clone2->Read().seq;


        auto v_shm_it1 = clone1->VSHMs().cbegin();
        auto v_shm_it2 = clone2->VSHMs().cbegin();

        size_t read1_length = seqan::length(read1);
        size_t read2_length = seqan::length(read2);

        TraverseReads(read1, read2,
                      0, 0,
                      clone1->CDR3Range().start_pos, clone2->CDR3Range().start_pos,
                      v_shm_it1, v_shm_it2,
                      clone1->VSHMs().cend(), clone2->VSHMs().cend(),
                      res_string);
        for (size_t i = clone1->CDR3Range().start_pos; i <= clone1->CDR3Range().end_pos; ++i) {
            res_string.push_back(read1[i]);
        }
        auto j_shm_it1 = clone1->JSHMs().cbegin();
        while (j_shm_it1 != clone1->JSHMs().cend() && j_shm_it1->read_nucl_pos <= clone1->CDR3Range().end_pos) {
            j_shm_it1++;
        }
        auto j_shm_it2 = clone2->JSHMs().cbegin();
        while (j_shm_it2 != clone2->JSHMs().cend() && j_shm_it2->read_nucl_pos <= clone2->CDR3Range().end_pos) {
            j_shm_it2++;
        }
        TraverseReads(read1, read2,
                      clone1->CDR3Range().end_pos + 1, clone2->CDR3Range().end_pos + 1,
                      read1_length, read2_length,
                      j_shm_it1, j_shm_it2,
                      clone1->JSHMs().cend(), clone2->JSHMs().cend(),
                      res_string);
        for (auto c : res_string) {
            seqan::append(res_dna5string, c);
        }
        core::Read res_read(read_name, res_dna5string, id);
        return res_read;
    }

    static void ParentReadReconstructor::TraverseReads(const seqan::String& read1,
                                                        const seqan::String& read2,
                                                        size_t read1_start_pos,
                                                        size_t read2_start_pos,
                                                        size_t read1_end_pos,
                                                        size_t read2_end_pos,
                                                        annotation_utils::GeneSegmentSHMs::SHMConstIterator shm_it1,
                                                        annotation_utils::GeneSegmentSHMs::SHMConstIterator shm_it2,
                                                        annotation_utils::GeneSegmentSHMs::SHMConstIterator shm_end1,
                                                        annotation_utils::GeneSegmentSHMs::SHMConstIterator shm_end2,
                                                        std::string& res_string) const {

        size_t read1_nucl_index = read1_start_pos;
        size_t read2_nucl_index = read2_start_pos;

        while (shm_it1 != shm_end1 &&
               shm_it2 != shm_end2 &&
               read1_nucl_index < read1_end_pos &&
               read2_nucl_index < read2_end_pos) {
            //todo: add cdr ranges assigning
            if (read1_nucl_index != shm_it1->read_nucl_pos && read2_nucl_index != shm_it2->read_nucl_pos) {
                VERIFY_MSG(read1[read1_nucl_index] == read2[read2_nucl_index],
                           "intersection_parent_constructor: unequal nucleotides on unmutated pos");
                res_string.push_back(read1[read1_nucl_index]);
                ++read1_nucl_index;
                ++read2_nucl_index;
                continue;
            }
            if (read1_nucl_index != shm_it1->read_nucl_pos) {
                if (shm_it2->shm_type == annotation_utils::SHMType::InsertionSHM) {
                    shm_it2++;
                    ++read2_nucl_index;
                    continue;
                }
                if (shm_it2->shm_type == annotation_utils::SHMType::DeletionSHM) {
                    shm_it2++;
                    continue;
                }
                if (shm_it2->shm_type == annotation_utils::SHMType::SubstitutionSHM) {
                    shm_it2++;
                    res_string.push_back(read1[read1_nucl_index]);
                    ++read1_nucl_index;
                    ++read2_nucl_index;
                    continue;
                }
                VERIFY_MSG(false, "intersection_parent_constructor: got unknown SHM type");
            }
            if (read2_nucl_index != shm_it2->read_nucl_pos) {
                if (shm_it1->shm_type == annotation_utils::SHMType::InsertionSHM) {
                    shm_it1++;
                    ++read1_nucl_index;
                    continue;
                }
                if (shm_it1->shm_type == annotation_utils::SHMType::DeletionSHM) {
                    shm_it1++;
                    continue;
                }
                if (shm_it1->shm_type == annotation_utils::SHMType::SubstitutionSHM) {
                    shm_it1++;
                    res_string.push_back(read2[read2_nucl_index]);
                    ++read1_nucl_index;
                    ++read2_nucl_index;
                    continue;
                }
                VERIFY_MSG(false, "intersection_parent_constructor: got unknown SHM type");
                continue;
            }
            if (shm_it1->shm_type == annotation_utils::SHMType::DeletionSHM &&
                shm_it2->shm_type == annotation_utils::SHMType::DeletionSHM) {
                shm_it1++;
                shm_it2++;
                auto new_shm = *shm_it1;
                new_shm.read_nucl_pos = res_string.size();

                continue;
            }
            if (shm_it1->shm_type == annotation_utils::SHMType::InsertionSHM &&
                shm_it2->shm_type == annotation_utils::SHMType::InsertionSHM) {
                if (read1[read1_nucl_index] == read2[read2_nucl_index]) { // identical insertion
                    res_string.push_back(read1[read1_nucl_index]);
                }
                shm_it1++;
                shm_it2++;
                ++read1_nucl_index;
                ++read2_nucl_index;
                continue;
            }
            if (shm_it1->shm_type == annotation_utils::SHMType::SubstitutionSHM &&
                shm_it2->shm_type == annotation_utils::SHMType::SubstitutionSHM) {
                if (read1[read1_nucl_index] == read2[read2_nucl_index]) { // identical substitution
                    res_string.push_back(read1[read1_nucl_index]);
                }
                else {
                    res_string.push_back(shm_it1->gene_nucl);
                }
                shm_it1++;
                shm_it2++;
                ++read1_nucl_index;
                ++read2_nucl_index;
                continue;
            }
            if (shm_it1->shm_type == annotation_utils::SHMType::InsertionSHM) {
                shm_it1++;
                ++read1_nucl_index;
                continue;
            }
            if (shm_it2->shm_type == annotation_utils::SHMType::InsertionSHM) {
                shm_it2++;
                ++read2_nucl_index;
                continue;
            }
            if (shm_it1->shm_type == annotation_utils::SHMType::SubstitutionSHM &&
                shm_it2->shm_type == annotation_utils::SHMType::DeletionSHM) {
                res_string.push_back(shm_it1->gene_nucl);
                shm_it1++;
                shm_it2++;
                ++read1_nucl_index;
                continue;
            }
            if (shm_it1->shm_type == annotation_utils::SHMType::DeletionSHM &&
                shm_it2->shm_type == annotation_utils::SHMType::SubstitutionSHM) {
                res_string.push_back(shm_it1->gene_nucl);
                shm_it1++;
                shm_it2++;
                ++read2_nucl_index;
                continue;
            }
        }
    }
}