#include "droplet_barcode.hpp"

bool operator<(const DropletBarcode& l, const DropletBarcode& r) {
    if(l.StrId() != r.StrId())
        return l.StrId() < r.StrId();
    return l.IntId() < r.IntId();
}

std::ostream& operator<<(std::ostream& out, const DropletBarcode& db) {
    out << db.StrId() << db.IntId();
    return out;
}

size_t DropletBarcodeHasher::operator()(const DropletBarcode &obj) const {
    return std::hash<std::string>()(obj.StrId()) * std::hash<unsigned long>()(obj.IntId());
}