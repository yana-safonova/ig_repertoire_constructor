#pragma once

#include <string>
#include <vector>
#include "../ig_tools/utils/string_tools.hpp"

template <typename ColumnTypes>
class ReportColumns {
public:
    ReportColumns() = default;

    ReportColumns(const std::vector<typename ColumnTypes::ColumnTypeEnum>& columns_, const std::string& csv_header_) :
            csv_header(csv_header_), columns(columns_) {}

    static ReportColumns CreateFromString(const std::string& columns_raw) {
        std::vector<std::string> col_strings = split(columns_raw, ',');
        std::vector<typename ColumnTypes::ColumnTypeEnum> columns;
        std::stringstream header;
        bool first = true;
        for (auto& s : col_strings) {
            boost::algorithm::trim(s);
            VERIFY_MSG(ColumnTypes::string_to_column_type.find(s) != ColumnTypes::string_to_column_type.end(), s);
            columns.push_back(ColumnTypes::string_to_column_type.at(s));
            if (!first) {
                header << "\t";
            }
            first = false;
            header << s;
        }
        return ReportColumns(columns, header.str());
    }

    std::string GetCsvHeader() const { return csv_header; }

    std::vector<typename ColumnTypes::ColumnTypeEnum> GetColumns() const { return columns; }

private:
    std::string csv_header;

    std::vector<typename ColumnTypes::ColumnTypeEnum> columns;
};
