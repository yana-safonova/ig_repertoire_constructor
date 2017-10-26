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
        VERIFY_MSG(false, "Could not find a supported " << type_name << " with name " << query << ". Supported " << type_name << "s are " << ss);
    }

    template <typename Context, typename Out>
    struct Column {
        const std::string name;
        const std::function<void(Out&, const Context&)> print;

        static Column<Context, Out> GetColumn(const std::string& name) {
            return find_element_by_name_or_report_failure<Column>(COLUMN_TYPES, name, "column");
        }

        static const std::vector<Column<Context, Out>> COLUMN_TYPES;
    };

    template <typename Context, typename Out>
    struct ColumnSet {
        const std::string name;
        const std::vector<Column<Context, Out>> columns;

        ColumnSet(const std::string& name, const std::vector<Column<Context, Out>>& columns) : name(name), columns(columns) {}

        void PrintCsvHeader(Out& out) const {
            PrintInfo(out, [&](Column<Context, Out> column) { out << column.name; });
        }

        void print(Out& out, const Context& context) const {
            PrintInfo(out, [&](Column<Context, Out> column) { column.print(out, context); });
        }

        static ColumnSet<Context, Out> GetPreset(const std::string& name) {
            return find_element_by_name_or_report_failure<ColumnSet>(PRESETS, name, "column set");
        }

        static ColumnSet<Context, Out> ParseColumns(const std::string& line) {
            std::vector<std::string> col_strings = split(line, ',');
            std::vector<Column<Context, Out>> columns;
            for (auto& name : col_strings) {
                boost::algorithm::trim(name);
                columns.push_back(Column<Context, Out>::GetColumn(name));

            }
            return ColumnSet{"custom", columns};
        }

        static const std::vector<ColumnSet<Context, Out>> PRESETS;

    private:
        void PrintInfo(Out& out, std::function<void(Column<Context, Out>)> print_info) const {
            bool first = true;
            for (const auto& column : columns) {
                if (!first) {
                    out << "\t";
                }
                first = false;
                print_info(column);
            }
            out << "\n";
        }
    };
};

