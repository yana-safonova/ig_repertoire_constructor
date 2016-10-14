#include <verify.hpp>
#include <path_helper.hpp>
#include <logger/logger.hpp>

#include "read_archive.hpp"

#include <seqan/seq_io.h>
#include <vector>


namespace core {
    std::ostream &operator<<(std::ostream &out, const Read &read) {
        out << "Name: " << read.name << ". Seq (" << seqan::length(read.seq) << "): " << read.seq;
        return out;
    }

//-----------------------------------------------------------------

    ReadArchive::ReadArchive(std::string fastq_file_fname) {
        ExtractFromFile(fastq_file_fname);
    }

    void ReadArchive::ExtractFromFile(std::string fastq_file_fname) {
        path::CheckFileExistenceFATAL(fastq_file_fname);
        std::vector <seqan::CharString> read_headers;
        std::vector <seqan::Dna5String> read_seqs;
        seqan::SeqFileIn seqFileIn_reads(fastq_file_fname.c_str());
        seqan::readRecords(read_headers, read_seqs, seqFileIn_reads);
        for (size_t i = 0; i < read_seqs.size(); i++) {
            std::string header = std::string(seqan::toCString(read_headers[i]));
            Read read(header, read_seqs[i], i);
            reads_.push_back(read);
            name_index_map_[header] = i;
        }
        INFO(size() << " reads were extracted from " << fastq_file_fname);
    }

    size_t ReadArchive::size() const {
        return reads_.size();
    }

    const Read& ReadArchive::operator[](size_t index) const {
        VERIFY_MSG(index < reads_.size(), "Index " << index << " exceeds archive size");
        return reads_[index];
    }

    size_t ReadArchive::GetIndexByReadName(std::string read_name) const {
        VERIFY_MSG(name_index_map_.find(read_name) != name_index_map_.end(),
                   "Read archive does not contain " << read_name);
        return name_index_map_.at(read_name);
    }

    const Read& ReadArchive::GetReadByName(std::string read_name) const {
        size_t index = GetIndexByReadName(read_name);
        return reads_[index];
    }

    void ReadArchive::FixSpacesInHeaders() {
        for(auto it = reads_.begin(); it != reads_.end(); it++)
            std::replace(it->name.begin(), it->name.end(), ' ', '_');
    }

    void ReadArchive::UpdateReadByIndex(size_t index, seqan::Dna5String new_seq) {
        VERIFY_MSG(index < reads_.size(), "Index " << index << " exceeds archive size");
        reads_[index].seq = new_seq;
    }
}
