#include <verify.hpp>

#include "germline_db_labeler.hpp"

namespace cdr_labeler {
    DbCDRLabeling GermlineDbLabeler::ComputeLabeling() {
        VERIFY_MSG(false, "Implement me");
        return DbCDRLabeling(gene_db_);
    }
}