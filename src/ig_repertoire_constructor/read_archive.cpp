#include "read_archive.hpp"
#include "verify.hpp"
#include "logger/logger.hpp"

namespace ig_repertoire_constructor {

io::SingleRead& ReadArchive::operator[](size_t number) {
    return reads_[number];
}

const io::SingleRead& ReadArchive::operator[](size_t number) const {
    return reads_[number];
}

size_t ReadArchive::GetReadNumberByReadName(const std::string & read_name) const {
    auto it = read_number_by_name_.find(read_name);
    VERIFY_MSG(it != read_number_by_name_.end(), "Read with name " << read_name << " not found");
    return it->second;
}

size_t ReadArchive::size() const {
    return reads_.size();
}

ReadArchive::iterator ReadArchive::begin() {
    return reads_.begin();
}

ReadArchive::iterator ReadArchive::end() {
    return reads_.end();
}

ReadArchive::const_iterator ReadArchive::begin() const {
    return reads_.begin();
}

ReadArchive::const_iterator ReadArchive::end() const {
    return reads_.end();
}

void ReadArchive::AddNewRead(const io::SingleRead & read) {
    VERIFY_MSG(read_number_by_name_.find(read.name()) == read_number_by_name_.end(), "Read with name " << read.name() << " already added");

    size_t read_number = reads_.size();
    reads_.push_back(read);
    read_number_by_name_.insert(make_pair(read.name(), read_number));
}

}
