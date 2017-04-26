//
// Created by Andrew Bzikadze on 4/12/17.
//

#pragma once

#include <memory>
#include "ig_simulator_utils.hpp"
#include "clonal_trees/tree/node.hpp"
#include <random>
#include "ig_simulator_config.hpp"

namespace ig_simulator {

class AbstractShmCreator {
protected:
    const size_t fix_left;
    const size_t fix_right;

public:
    AbstractShmCreator() = delete;
    AbstractShmCreator(const AbstractShmCreator&) = delete;
    AbstractShmCreator(AbstractShmCreator&&) = delete;
    AbstractShmCreator& operator=(const AbstractShmCreator&) = delete;
    AbstractShmCreator& operator=(AbstractShmCreator&&) = delete;

    explicit AbstractShmCreator(const vj_finder::VJFinderConfig& config):
        fix_left(config.algorithm_params.fix_crop_fill_params.fix_left),
        fix_right(config.algorithm_params.fix_crop_fill_params.fix_right)
    { }

    virtual ~AbstractShmCreator() { }

    virtual Node::SHM_Vector GenerateSHM_Vector(const std::string&) const = 0;
};

using AbstractShmCreatorCPtr = std::unique_ptr<AbstractShmCreator>;


class PoissonShmCreator final : public AbstractShmCreator {
private:
    mutable std::poisson_distribution<size_t> distribution;

public:
    PoissonShmCreator(const vj_finder::VJFinderConfig& vjf_config,
                      double lambda):
        AbstractShmCreator(vjf_config),
        distribution(check_numeric_positive(lambda))
    { }

    PoissonShmCreator(const vj_finder::VJFinderConfig& vjf_config,
                      const SHM_CreatorParams::PoissonCreatorParams& config):
        PoissonShmCreator(vjf_config, config.lambda)
    { }


    Node::SHM_Vector GenerateSHM_Vector(const std::string&) const override;
};

AbstractShmCreatorCPtr get_shm_creator(const vj_finder::VJFinderConfig&, const SHM_CreatorParams&);

} // End namespace ig_simulator