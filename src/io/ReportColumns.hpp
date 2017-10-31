#pragma once

#include <string>
#include <vector>
#include "../ig_tools/utils/string_tools.hpp"

namespace ReportColumns {
    template <typename T, typename TList>
    T find_element_by_name_or_report_failure(const TList& list, const std::string& query, const std::string& type_name) {
        std::stringstream ss;
        bool first = true;
        for (const auto& element : list) {
            if (element.name == query) {
                return element;
            }
            if (!first) ss << ", ";
            first = false;
            ss << element.name;
        }
        VERIFY_MSG(false, "Could not find a supported " << type_name << " with name " << query << ". Supported " << type_name << "s are " << ss.str());
    }

    template <typename Context>
    struct Column {
        const std::string name;
        const std::function<void(std::basic_ostream<char>&, const Context&)> print;

        static Column<Context> GetColumn(const std::string& name) {
            return find_element_by_name_or_report_failure<Column>(COLUMN_TYPES, name, "column");
        }

        static const std::vector<Column<Context>> COLUMN_TYPES;
    };

    template <typename Context>
    struct ColumnSet {
        const std::string name;
        const std::vector<Column<Context>> columns;

        ColumnSet(const std::string& name, const std::vector<Column<Context>>& columns) : name(name), columns(columns) {}

        void PrintCsvHeader(std::basic_ostream<char>& out, const std::string& separator = "\t") const {
            PrintInfo(out, [&](Column<Context> column) { out << column.name; }, separator);
        }

        void Print(std::basic_ostream<char>& out, const Context& context) const {
            PrintInfo(out, [&](Column<Context> column) { column.print(out, context); });
        }

        static ColumnSet<Context> GetPreset(const std::string& name) {
            return find_element_by_name_or_report_failure<ColumnSet>(PRESETS, name, "column set");
        }

        static ColumnSet<Context> ParseColumns(const std::string& line) {
            std::vector<std::string> col_strings = split(line, ',');
            std::vector<Column<Context>> columns;
            for (auto& name : col_strings) {
                boost::algorithm::trim(name);
                columns.push_back(Column<Context>::GetColumn(name));

            }
            return ColumnSet{"custom", columns};
        }

        static ColumnSet<Context> ChooseColumns(const std::string& preset_name, const std::string& columns) {
            if (!columns.empty()) {
                return ParseColumns(columns);
            }
            return GetPreset(preset_name);
        };

        static const std::vector<ColumnSet<Context>> PRESETS;

    private:
        void PrintInfo(std::basic_ostream<char>& out, const std::function<void(Column<Context>)>& print_info, const std::string& separator = "\t") const {
            bool first = true;
            for (const auto& column : columns) {
                if (!first) {
                    out << separator;
                }
                first = false;
                print_info(column);
            }
            out << "\n";
        }
    };
};

