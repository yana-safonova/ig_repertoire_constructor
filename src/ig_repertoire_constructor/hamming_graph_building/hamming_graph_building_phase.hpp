#pragma once

#include "stage.hpp"
#include "read_archive.hpp"
#include "ig_repertoire_constructor.hpp"
#include <fstream>

#include "hamming_graph_builder.hpp"

namespace ig_repertoire_constructor {

class HammingGraphBuildingPhase : public IgRepertoireConstructor::Phase {
public:
    HammingGraphBuildingPhase() : IgRepertoireConstructor::Phase("Building Hamming graph on kmers", "hamming_graph_building") { }

    void run(const char*) {
        auto read_archive_ptr = LoadReadArchive();
        storage().SetReadArchive(read_archive_ptr);

        auto components_ptr = HammingGraphBuilder().Build(*read_archive_ptr);
        storage().SetHammingComponents(components_ptr);
    }

    void load(const std::string &load_from,
            const char* prefix) {
        std::string file_name = path::append_path(load_from, prefix) + ".hamming_components";
        INFO("Loading current state from " << file_name);

        auto read_archive_ptr = LoadReadArchive();
        storage().SetReadArchive(read_archive_ptr);

        MMappedReader ifs(file_name.c_str(), std::ios_base::in);
        VERIFY(ifs.good());

        auto components_ptr = std::make_shared <std::vector <std::vector <size_t> > >();
        size_t count;
        ifs.read((char*)&count, sizeof(count));
        components_ptr->resize(count);
        for (auto & component : *components_ptr) {
            size_t size;
            ifs.read((char*)&size, sizeof(size));
            component.resize(size);
            ifs.read((char*)component.data(), sizeof(size_t) * size);
        }

        INFO(components_ptr->size() << " Hamming components have been loaded");
        storage().SetHammingComponents(components_ptr);
    }

    void save(const std::string & save_to,
            const char* prefix) const {
        std::string file_name = path::append_path(save_to, prefix) + ".hamming_components";
        INFO("Saving current state to " << file_name);

        auto read_archive_ptr = storage().GetReadArchivePtr();
        auto components_ptr = storage().GetHammingComponentsPtr();

        // TODO save ReadArchive?

        std::ofstream ofs(file_name.c_str(), std::ios_base::out);
        VERIFY(ofs.good());

        size_t count = components_ptr->size();
        ofs.write((char*)&count, sizeof(count));
        for (const auto & component : *components_ptr) {
            size_t size = component.size();
            ofs.write((char*)&size, sizeof(size));
            ofs.write((char*)component.data(), sizeof(size_t) * size);
        }

        VERIFY(!ofs.fail());
        ofs.close();
        INFO(components_ptr->size() << " Hamming components have been saved");
    }

    virtual ~HammingGraphBuildingPhase() { }
};

}
