#include "vj_alignment_info.hpp"
#include <boost/algorithm/string.hpp>
#include <../ig_tools/utils/string_tools.hpp>

using namespace std;
using namespace boost;

void VJAlignmentInfo::AddReadName(const vector<string> &splits) {
    string read_name = splits[VJFinderMagicConsts::read_name_index];
    read_names_.push_back(read_name);
}

void VJAlignmentInfo::AddVAlignment(const std::vector<std::string> &splits) {
    size_t start_read_pos = string_to_number<size_t>(splits[VJFinderMagicConsts::first_v_match_pos_on_read]);
    size_t end_read_pos = string_to_number<size_t>(splits[VJFinderMagicConsts::last_v_match_pos_on_read]);
    size_t start_gene_pos = string_to_number<size_t>(splits[VJFinderMagicConsts::first_v_match_pos_on_gene]);
    size_t end_gene_pos = string_to_number<size_t>(splits[VJFinderMagicConsts::last_v_match_pos_on_gene]);
    string v_gene_name = splits[VJFinderMagicConsts::v_gene_name_index];
    v_segments_.push_back(IgGeneAlignmentPositions(
                    AlignmentPositions(pair<size_t, size_t>(start_read_pos, end_read_pos),
                              pair<size_t, size_t>(start_gene_pos, end_gene_pos)),
                    v_gene_db_.GetByName(v_gene_name)));
}

void VJAlignmentInfo::AddJAlignment(const std::vector<std::string> &splits) {
    size_t start_read_pos = string_to_number<size_t>(splits[VJFinderMagicConsts::first_j_match_pos_on_read]);
    size_t end_read_pos = string_to_number<size_t>(splits[VJFinderMagicConsts::last_j_match_pos_on_read]);
    size_t start_gene_pos = string_to_number<size_t>(splits[VJFinderMagicConsts::first_j_match_pos_on_gene]);
    size_t end_gene_pos = string_to_number<size_t>(splits[VJFinderMagicConsts::last_j_match_pos_on_gene]);
    string j_gene_name = splits[VJFinderMagicConsts::j_gene_name_index];
    cout << j_gene_name << endl;
    j_segments_.push_back(IgGeneAlignmentPositions(
            AlignmentPositions(pair<size_t, size_t>(start_read_pos, end_read_pos),
                      pair<size_t, size_t>(start_gene_pos, end_gene_pos)),
            j_gene_db_.GetByName(j_gene_name)));
}

void VJAlignmentInfo::ParseLine(std::string line) {
    vector<std::string> splits;
    split(splits, line, is_any_of(" ,"), token_compress_on);
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

IgGeneAlignmentPositions VJAlignmentInfo::GetVAlignmentByIndex(size_t index) const {
    assert(index < v_segments_.size());
    return v_segments_[index];
}

IgGeneAlignmentPositions VJAlignmentInfo::GetJAlignmentByIndex(size_t index) const {
    assert(index < j_segments_.size());
    return j_segments_[index];
}

std::string VJAlignmentInfo::GetReadNameByIndex(size_t index) const {
    assert(index < read_names_.size());
    return read_names_[index];
}

std::ostream& operator<<(std::ostream& out, const VJAlignmentInfo& obj) {
    for(size_t i = 0; i < obj.size(); i++) {
        out << "Read name: " << obj.GetReadNameByIndex(i) << endl;
        out << "V. " << obj.GetVAlignmentByIndex(i) << endl;
        out << "J. " << obj.GetJAlignmentByIndex(i) << endl;
        out << "----" << endl;
    }
    return out;
}