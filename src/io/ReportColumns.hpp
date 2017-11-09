#pragma once

#include <string>
#include <vector>
#include "../ig_tools/utils/string_tools.hpp"

namespace ReportColumns {
    template <typename T, typename TList>
    boost::optional<T> get_element_by_name(const TList& list, const std::string& query) {
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
        return boost::optional<T>();
    }

    template <typename Context>
    struct Column {
        const std::string name;
        const std::function<void(std::ostream&, const Context&)> print;

        static boost::optional<Column<Context>> GetColumn(const std::string& name) {
            return get_element_by_name<Column>(COLUMN_TYPES, name);
        }

        static const std::vector<Column<Context>> COLUMN_TYPES;
    };

    template <typename Context>
    struct ColumnSet {
        const std::string name;
        const std::vector<Column<Context>> columns;

        ColumnSet(const std::string& name, const std::vector<Column<Context>>& columns) : name(name), columns(columns) {}

        void PrintCsvHeader(std::ostream& out, const std::string& separator = "\t") const {
            PrintInfo(out, [&](Column<Context> column) { out << column.name; }, separator);
        }

        void Print(std::ostream& out, const Context& context) const {
            PrintInfo(out, [&](Column<Context> column) { column.print(out, context); });
        }

        static boost::optional<ColumnSet<Context>> GetPreset(const std::string& name) {
            return get_element_by_name<ColumnSet>(PRESETS, name);
        }

        static boost::optional<ColumnSet<Context>> ParseColumns(const std::string& line) {
            std::vector<std::string> col_strings = split(line, ',');
            std::vector<Column<Context>> columns;
            for (auto& name : col_strings) {
                boost::algorithm::trim(name);
                const auto column = Column<Context>::GetColumn(name);
                if (!column) return boost::optional<ColumnSet<Context>>();
                columns.push_back(column.value());
            }
            return ColumnSet{"custom", columns};
        }

        static ColumnSet<Context> ChooseColumns(const std::string& preset_name, const std::string& columns) {
            if (!columns.empty()) {
                const auto columns_from_list = ParseColumns(columns);
                if (!columns_from_list) {
                    FATAL_ERROR("Could not parse column list from " << columns);
                }
                return columns_from_list.value();
            }
            const auto columns_from_preset = GetPreset(preset_name);
            if (!columns_from_preset) {
                FATAL_ERROR("Unknown preset name: " << preset_name);
            }
            return columns_from_preset.value();
        };

        static const std::vector<ColumnSet<Context>> PRESETS;

    private:
        void PrintInfo(std::ostream& out, const std::function<void(Column<Context>)>& print_info, const std::string& separator = "\t") const {
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

