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
}