#pragma once

#include <string>
#include <ostream>

class DropletBarcode {
    std::string str_id_;
    size_t int_id_;

public:
    DropletBarcode(std::string str_id, size_t int_id) :
            str_id_(str_id),
            int_id_(int_id) { }

    std::string StrId() const { return str_id_; }

    size_t IntId() const { return int_id_; }

    bool IsEmpty() const { return str_id_ == ""; }

    bool operator ==(const DropletBarcode &db) const { return str_id_ == db.str_id_ and int_id_ == db.int_id_; }
};

bool operator<(const DropletBarcode& l, const DropletBarcode& r);

std::ostream& operator<<(std::ostream& out, const DropletBarcode& db);

struct DropletBarcodeHasher {
    size_t operator()(const DropletBarcode &obj) const;
};