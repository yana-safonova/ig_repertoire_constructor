#include "shm_annotation.hpp"

namespace annotation_utils {
    std::ostream& operator<<(std::ostream& out, const SHMType &shm_type) {
        if(shm_type == SHMType::SubstitutionSHM)
            out << "S";
        else if(shm_type == SHMType::DeletionSHM)
            out << "D";
        else if(shm_type == SHMType::InsertionSHM)
            out << "I";
        else
            out << "Unknown_SHM";
        return out;
    }

    void SHM::ComputeType() {
        if(gene_nucl == '-')
            shm_type = SHMType::InsertionSHM;
        else if(read_nucl == '-')
            shm_type = SHMType::DeletionSHM;
        else
            shm_type = SHMType::SubstitutionSHM;
    }

    bool SHM::operator==(const SHM &shm) const {
        if(gene_nucl_pos != shm.gene_nucl_pos)
            return false;
        if(gene_nucl != shm.gene_nucl)
            return false;
        if(read_nucl != shm.read_nucl)
            return false;
//        if(read_aa != shm.read_aa)
//            return false;
        return shm_type == shm.shm_type;
    }

    bool SHM::operator!=(const SHM &shm) const {
        return !(*this == shm);
    }

    bool operator<(const SHM &left, const SHM &right) {
        if(left.gene_nucl_pos != right.gene_nucl_pos)
            return left.gene_nucl_pos < right.gene_nucl_pos;
        if(left.read_aa != right.read_aa)
            return left.read_aa < right.read_aa;
        //if(left.read_aa != right.read_aa)
        return left.read_nucl < right.read_nucl;
        //return left.gene_nucl < right.gene_nucl;
    }

    bool TrivialSHMComparator::operator()(const SHM &shm1, const SHM &shm2) {
        if(shm1.gene_nucl_pos != shm2.gene_nucl_pos)
            return shm1.gene_nucl_pos < shm2.gene_nucl_pos;
        if(shm1.read_nucl_pos != shm2.read_nucl_pos)
            return shm1.read_nucl_pos < shm2.read_nucl_pos;
        if(shm1.gene_nucl != shm2.gene_nucl)
            return shm1.gene_nucl < shm2.gene_nucl;
        if(shm1.read_nucl != shm2.read_nucl)
            return shm1.read_nucl < shm2.read_nucl;
        if(shm1.gene_aa != shm2.gene_aa)
            return shm1.gene_aa < shm2.gene_aa;
        return shm1.read_aa < shm2.read_aa;
    }

    std::ostream& operator<<(std::ostream &out, const SHM& shm) {
        out << shm.gene_nucl_pos << " - " << shm.read_nucl_pos << ", " << shm.gene_nucl << "->" << shm.read_nucl << ", " <<
        shm.gene_aa << "->" << shm.read_aa;
        return out;
    }

    void GeneSegmentSHMs::CheckConsistencyFatal(SHM shm) {
        if(size() == 0)
            return;
        SHM last_shm = shms_[size() - 1];
        VERIFY_MSG(!(last_shm.read_nucl_pos > shm.read_nucl_pos or last_shm.gene_nucl_pos > shm.gene_nucl_pos),
                   "Order of SHMs " << last_shm << " and " << shm << " are not consistent");
        VERIFY_MSG(last_shm.read_nucl_pos != shm.read_nucl_pos or last_shm.gene_nucl_pos != shm.gene_nucl_pos,
                   "SHMs " << last_shm << " and " << shm << " have identical read and gene positions");
    }

    void GeneSegmentSHMs::AddSHM(SHM shm) {
        CheckConsistencyFatal(shm);
        shms_.push_back(shm);
    }

    const SHM& GeneSegmentSHMs::operator[](size_t index) const {
        VERIFY_MSG(index < size(), "Index " << index << " exceeds number of SHMs");
        return shms_[index];
    }

    bool GeneSegmentSHMs::operator==(const SHM &) const {
        VERIFY_MSG(false, "Implement me!");
        return false;
    }

    std::ostream& operator<<(std::ostream &out, const GeneSegmentSHMs& shms) {
        for(auto it = shms.cbegin(); it != shms.cend(); it++)
            out << *it << std::endl;
        return out;
    }

    void SHM::AppendInMixcrFormat(std::ostream& out) const {
        out << shm_type;
        if (shm_type != InsertionSHM) {
            out << gene_nucl;
        }
        out << gene_nucl_pos;
        if (shm_type != DeletionSHM) {
            out << read_nucl;
        }
    }

    void GeneSegmentSHMs::AppendInMixcrFormat(std::ostream& out) const {
        for (auto shm_it = cbegin(); shm_it != cend(); ++ shm_it) {
            shm_it->AppendInMixcrFormat(out);
        }
    }
}