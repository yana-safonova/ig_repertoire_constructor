//
// Created by Andrew Bzikadze on 4/9/17.
//

#pragma once

#include "clonal_trees/tree/tree.hpp"
#include "base_repertoire/metaroot_cluster/metaroot_cluster.hpp"

namespace ig_simulator {

class Forest {
private:
    const MetarootCluster* metaroot_cluster;
    const std::vector<Tree> trees;

public:
    Forest(const MetarootCluster* const metaroot_cluster,
           std::vector<Tree>&& trees = {}) noexcept:
        metaroot_cluster(metaroot_cluster),
        trees(trees)
    { }

    Forest(const Forest&) = default;
    Forest(Forest&&) = default;

    Forest& operator=(const Forest&) = default;
    Forest& operator=(Forest&&) = default;

    const MetarootCluster* GetMetarootCluster() const { return metaroot_cluster; }
    const std::vector<Tree>& Trees() const { return trees; }

    size_t Size() const { return trees.size(); }

    friend std::ostream& operator<<(std::ostream&, const Forest&);
};

std::ostream& operator<<(std::ostream& out, const Forest&);

using ForestStorage = std::vector<Forest>;

} // End namespace ig_simulator
