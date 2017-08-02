#include "shm_comparator.hpp"

namespace annotation_utils {
    bool SHMComparator::SHMsAreEqual(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2) {
        if(shms1.size() != shms2.size())
            return false;
        for(size_t i = 0; i < shms1.size(); i++)
            if(shms1[i] != shms2[i])
                return false;
        return true;
    }

    static std::vector<SHM> get_following_insertions_block(GeneSegmentSHMs::SHMConstIterator it,
                                                           GeneSegmentSHMs::SHMConstIterator end) {
        std::vector<SHM> res;
        do {
            res.push_back(*it);
            ++it;
        } while (it != end &&
                 it->shm_type == SHMType::InsertionSHM &&
                 it->gene_nucl_pos == res.back().gene_nucl_pos &&
                 it->read_nucl_pos == res.back().read_nucl_pos + 1);
        return res;
    }
    static bool insertion_blocks_are_equal(const std::vector<SHM>& block1,
                                           const std::vector<SHM>& block2) {
        if (block1.size() != block2.size()) {
            return false;
        }
        for (size_t i = 0; i < block1.size(); ++i) {
            const auto& shm1 = block1[i];
            const auto& shm2 = block2[i];
            if (shm1.read_nucl != shm2.read_nucl) {
                return false;
            }
        }
        return true;
    }

    bool SHMComparator::SHMs1AreNestedInSHMs2(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2) {
        size_t index2 = 0;
        for(auto it1 = shms1.cbegin(); it1 != shms1.cend(); it1++) {
            if (it1->shm_type == SHMType::InsertionSHM) {
                continue;
            }
            bool shm_found = false;
            for(size_t i = index2; i < shms2.size(); i++)
                if(*it1 == shms2[i]) {
                    index2 = i + 1;
                    shm_found = true;
                    break;
                }
            if(!shm_found)
                return false;
        }
        return AllSHMs1InsertionBlocksArePresentedInSHMs2(shms1, shms2);
//        return SHMsInsertionBlocksAreEqual(shms1, shms2);
    }

    bool SHMComparator::AllSHMs1InsertionBlocksArePresentedInSHMs2(GeneSegmentSHMs shms1,
                                                                   GeneSegmentSHMs shms2) {
        size_t index2 = 0;
        for (auto it1 = shms1.cbegin(); it1 != shms1.cend(); it1++) {
            bool shm_found = false;
            if (it1->shm_type != SHMType::InsertionSHM) {
                continue;
            }
            auto shm1_block = get_following_insertions_block(it1, shms1.cend());
            for(size_t i = index2; i < shms2.size(); i++) {
                if (*it1 == shms2[i]) {
                    auto shm2_block = get_following_insertions_block(shms2.cbegin() + i, shms2.cend());
                    if (insertion_blocks_are_equal(shm1_block, shm2_block)) {
                        shm_found = true;
                        it1 += shm1_block.size() - 1;
                        index2 = i + shm2_block.size();
                        break;
                    }
                    return false;
                }
            }
            if(!shm_found)
                return false;
        }
        return true;
    }
    bool SHMComparator::SHMsInsertionBlocksAreEqual(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2)
    {
        return AllSHMs1InsertionBlocksArePresentedInSHMs2(shms1, shms2) &&
               AllSHMs1InsertionBlocksArePresentedInSHMs2(shms2, shms1);
    }

    size_t SHMComparator::GetNumberOfIntersections(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2) {
        size_t num_shared_shms = 0;
        size_t index2 = 0;
        for(auto it1 = shms1.cbegin(); it1 != shms1.cend(); it1++) {
            for(size_t i = index2; i < shms2.size(); i++)
                if(*it1 == shms2[i]) {
                    num_shared_shms++;
                    index2 = i + 1;
                    break;
                }
        }
        return num_shared_shms;
    }

    bool SHMComparator::AddedSHMsAreSynonimous(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2) {
        size_t index1 = 0;
        for(auto it2 = shms2.cbegin(); it2 != shms2.cend(); it2++) {
            bool shm_found = false;
            for(size_t i = index1; i < shms1.size(); i++)
                if(*it2 == shms1[i]) {
                    index1 = i + 1;
                    shm_found = true;
                    break;
                }
            if(!shm_found && !it2->IsSynonymous()) {
                return false;
            }
        }
        return true;
    }

    bool SHMComparator::AllAddedSHMs1HaveIdenticallyPositionedSHMs2(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2) {
        size_t index1 = 0;
        for(auto it2 = shms2.cbegin(); it2 != shms2.cend(); it2++) {
            bool shm_found = false;
            for(size_t i = index1; i < shms1.size(); i++)
                if(it2->gene_nucl_pos == shms1[i].gene_nucl_pos &&
                   it2->gene_nucl == shms1[i].gene_nucl &&
                   it2->shm_type == shms1[i].shm_type) {
                    index1 = i + 1;
                    shm_found = true;
                    break;
                }
            if(!shm_found) {
                return false;
            }
        }
        return true;
    }

    bool SHMComparator::IndividualSHMsAreIdenticallyPositioned(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2) {
        return SHMsInsertionBlocksAreEqual(shms1, shms2) &&
               AllAddedSHMs1HaveIdenticallyPositionedSHMs2(shms1, shms2);
    }

    // return a vector of SHMs that appeared in SHM2, but not presented in SHM1
    std::vector<SHM> SHMComparator::GetAddedSHMs(GeneSegmentSHMs shms1, GeneSegmentSHMs shms2) {
        std::vector<SHM> added_shms;
        for(auto it2 = shms2.cbegin(); it2 != shms2.cend(); it2++) {
            bool shm_found = false;
            for(auto it1 = shms1.cbegin(); it1 != shms1.cend(); it1++) {
                if(*it1 == *it2) {
                    shm_found = true;
                    break;
                }
            }
            if(!shm_found) {
                added_shms.push_back(*it2);
            }
        }
        return added_shms;
    }
}
