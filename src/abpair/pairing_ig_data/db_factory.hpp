#pragma once

#include <string>
#include <unordered_map>
#include "pairing_primitives/droplet_barcode.hpp"

class DbFactory {
    std::unordered_map<std::string, size_t> filename_index_map_;

public:
    DropletBarcode GetDropletBarcodeByFilename(std::string db_str, std::string filename);
};