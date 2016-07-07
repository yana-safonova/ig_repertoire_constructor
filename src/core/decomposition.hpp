#pragma once

#include <vector>
#include <set>
#include <fstream>
#include <memory>

namespace core{
    class Decomposition {
        // input params
        size_t num_vertices_;

        // decomposition fields
        std::vector<size_t> vertex_class_;
        std::vector<std::set<size_t> > decomposition_classes_;

        // number of all classes: real and removed
        size_t num_classes_;

        void InitializeVertexClasses();

        std::vector<size_t> ReadClassIdsFromIfstream(std::ifstream &in);

        void AddNewClass();

        bool ClassIsValid(size_t class_id) { return class_id != size_t(-1); }

        void RemoveVertex(size_t vertex);

    public:
        Decomposition(size_t num_vertices) :
                num_vertices_(num_vertices) {
            InitializeVertexClasses();
        }

        Decomposition(std::string decomposition_filename);

        void SetClass(size_t vertex, size_t class_id);

        void AddDecomposition(std::shared_ptr<Decomposition> decomposition);

        const std::set<size_t>& GetClass(size_t index) const;

        const std::set<size_t>& LastClass() const { return GetClass(Size() - 1); }

        size_t ClassSize(size_t index) const { return GetClass(index).size(); }

        size_t LastClassSize() const { return LastClass().size(); }

        size_t LastClassId() const { return Size() - 1; }

        size_t VertexNumber() const { return vertex_class_.size(); }

        size_t GetVertexClass(size_t vertex) const;

        bool VertexClassIsInitialized(size_t vertex) {
            return ClassIsValid(GetVertexClass(vertex));
        }

        size_t Size() const { return decomposition_classes_.size(); }

        void SaveTo(std::string output_fname);

        bool LastClassContains(size_t vertex) const {
            return LastClass().find(vertex) != LastClass().end();
        }

        size_t MaxClassSize();

        bool IsTrivial() { return Size() == 1; }

        size_t NextClassId() { return num_classes_; }

    //private:
    //    DECL_LOGGER("Decomposition");
    };

    typedef std::set<size_t> DecompositionClass;

    std::ostream& operator<<(std::ostream &out, const Decomposition &hg_decomposition);

    typedef std::shared_ptr<Decomposition> DecompositionPtr;
}

