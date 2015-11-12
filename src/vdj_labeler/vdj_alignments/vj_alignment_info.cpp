#include "vj_alignment_info.hpp"
#include "../../ig_tools/utils/string_tools.hpp"

using namespace std;

void VJAlignmentInfo::AddReadName(const vector<string> &splits) {
    string read_name = splits[0].substr(0, splits[0].size() - 1);
    read_names_.push_back(read_name);
}

void VJAlignmentInfo::AddVAlignment(const std::vector<std::string> &splits) {
    size_t start_query_pos = string_to_number<size_t>(splits[1].substr(0, splits[1].size() - 1));
    size_t end_query_pos = string_to_number<size_t>(splits[2].substr(0, splits[2].size() - 1));
    string v_gene_name = splits[4].substr(0, splits[4].size() - 1);
    v_segments_.push_back(
            IgGeneAlignment(
                    Alignment(pair<size_t, size_t>(start_query_pos, end_query_pos),
                              pair<size_t, size_t>(0, 0)),
                    v_gene_db_.GetByName(v_gene_name)));
}

void VJAlignmentInfo::AddJAlignment(const std::vector<std::string> &splits) {
    size_t start_query_pos = string_to_number<size_t>(splits[5].substr(0, splits[5].size() - 1));
    size_t end_query_pos = string_to_number<size_t>(splits[6].substr(0, splits[6].size() - 1));
    string j_gene_name = splits[8].substr(0, splits[8].size() - 1);
    j_segments_.push_back(
            IgGeneAlignment(
                    Alignment(pair<size_t, size_t>(start_query_pos, end_query_pos),
                              pair<size_t, size_t>(0, 0)),
                    j_gene_db_.GetByName(j_gene_name)));
}

void VJAlignmentInfo::ParseLine(std::string line) {
    vector<std::string> splits = split(line, ' ');
    AddReadName(splits);
    AddVAlignment(splits);
    AddJAlignment(splits);
}

void VJAlignmentInfo::ExtractAlignment(std::string filename) {
    ifstream input(filename);
    if(!input.good()) {
        cout << "CSV file " << filename << "containing VJ labeling was not found";
        return;
    }
    size_t num_row = 0;
    while(!input.eof()) {
        num_row++;
        std::string tmp;
        getline(input, tmp);
        if(tmp == "")
            break;
        if(num_row == 1)
            continue;
        ParseLine(tmp);
    }
}

IgGeneAlignment VJAlignmentInfo::GetVAlignmentByIndex(size_t index) const {
    assert(index < v_segments_.size());
    return v_segments_[index];
}

IgGeneAlignment VJAlignmentInfo::GetJAlignmentByIndex(size_t index) const {
    assert(index < j_segments_.size());
    return j_segments_[index];
}

std::string VJAlignmentInfo::GetReadNameByIndex(size_t index) const {
    assert(index < read_names_.size());
    return read_names_[index];
}

std::ostream& operator<<(std::ostream& out, VJAlignmentInfo& obj) {
    for(size_t i = 0; i < obj.size(); i++) {
        out << obj.GetReadNameByIndex(i) << " " << obj.GetVAlignmentByIndex(i).alignment.query_pos.second << " " <<
                obj.GetJAlignmentByIndex(i).alignment.query_pos.first << endl;
    }
    return out;
}