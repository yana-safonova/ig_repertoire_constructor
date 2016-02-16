#include "db_factory.hpp"

DropletBarcode DbFactory::GetDropletBarcodeByFilename(std::string db_str, std::string filename) {
    size_t fname_index = size_t(-1);
    if(filename_index_map_.find(filename) == filename_index_map_.end()) {
        size_t old_map_size = filename_index_map_.size();
        filename_index_map_[filename] = old_map_size;
        fname_index = old_map_size;
    }
    else
        fname_index = filename_index_map_[filename];
    return DropletBarcode(db_str, fname_index);
}