//
// Created by Andrew Bzikadze on 3/15/17.
//

#include <chrono>

#include <clonal_trees/tree_creator/tree_creator.hpp>
#include <clonal_trees/tree_creator/exporters.hpp>
#include "ig_simulator_launch.hpp"
#include "base_repertoire/base_repertoire_simulator.hpp"
#include "clonal_trees/tree_creator/forest_storage_creator.hpp"

using namespace germline_utils;

namespace ig_simulator {

germline_utils::ChainType IgSimulatorLaunch::GetLaunchChainType() const {
    auto v_chain_type = germline_utils::LociParam::ConvertIntoChainTypes(config_.germline_params.loci);
    VERIFY_MSG(v_chain_type.size() == 1, "Only specific chain type is allowed");
    return v_chain_type[0];
}

std::vector<germline_utils::CustomGeneDatabase>
IgSimulatorLaunch::GetDB(const germline_utils::ChainType chain_type) const
{
    GermlineDbGenerator db_generator(config_.io_params.input_params.germline_input,
                                     config_.germline_params);
    INFO("Generation of DB for variable segments...");
    germline_utils::CustomGeneDatabase v_db = db_generator.GenerateVariableDb();
    INFO("Generation of DB for diversity segments...");
    germline_utils::CustomGeneDatabase d_db = db_generator.GenerateDiversityDb();
    INFO("Generation of DB for join segments...");
    germline_utils::CustomGeneDatabase j_db = db_generator.GenerateJoinDb();

    std::vector<germline_utils::CustomGeneDatabase> db;
    db.emplace_back(std::move(v_db));
    if (chain_type.IsVDJ())
        db.emplace_back(std::move(d_db));
    db.emplace_back(std::move(j_db));
    return db;
}

BaseRepertoire
IgSimulatorLaunch::GetBaseRepertoire(const germline_utils::ChainType chain_type,
                                     std::vector<germline_utils::CustomGeneDatabase>& db) const
{
    INFO("== Base Repertoire starts ==");
    BaseRepertoireSimulator base_repertoire_simulator{config_.simulation_params.base_repertoire_params,
                                                      chain_type,
                                                      db};
    auto base_repertoire =
        base_repertoire_simulator.Simulate(config_.simulation_params.base_repertoire_params.number_of_metaroots);
    std::ofstream base_repertoire_fasta;
    std::ofstream base_repertoire_info;
    base_repertoire_fasta.open(path::append_path(config_.io_params.output_params.output_dir,
                                                 config_.io_params.output_params.base_repertoire_filename));
    base_repertoire_info.open(path::append_path(config_.io_params.output_params.output_dir,
                                                config_.io_params.output_params.base_repertoire_info));
    print_base_repertoire(base_repertoire, base_repertoire_fasta, base_repertoire_info);
    base_repertoire_fasta.close();
    base_repertoire_info.close();
    INFO("== Base Repertoire ends ==");
    return base_repertoire;
}

template<class PoolManager>
ForestStorage IgSimulatorLaunch::__GetForestStorage(const BaseRepertoire& base_repertoire) const
{
    INFO("== Forest Storage generation starts ==");
    const auto& vjf_config = config_.simulation_params.base_repertoire_params.metaroot_simulation_params.
                             cdr_labeler_config.vj_finder_config;
    ForestStorageCreator forest_storage_creator(vjf_config,
                                                config_.simulation_params.clonal_tree_simulator_params);
    auto forest_storage = forest_storage_creator.GenerateForest<PoolManager>(base_repertoire);
    INFO("== Forest Storage generation ends ==");

    INFO("== Forest Storage export starts ==");
    INFO("== Full and filtered pool export start");
    std::ofstream full, included;
    full.open(path::append_path(config_.io_params.output_params.output_dir,
                                config_.io_params.output_params.full_pool));
    included.open(path::append_path(config_.io_params.output_params.output_dir,
                                    config_.io_params.output_params.filtered_pool));
    ForestStorageExporter(forest_storage, full, included);
    full.close();
    included.close();
    INFO("== Full and filtered pool export ends");

    INFO("== Edge lists export starts");
    EdgeListsExporters(forest_storage, config_.io_params.output_params);
    INFO("== Edge lists export ends");
    INFO("== Forest Storage export ends ==");
    return forest_storage;
}

ForestStorage IgSimulatorLaunch::GetForestStorage(const BaseRepertoire& base_repertoire) const
{
    const auto& pool_manager_strategy = config_.simulation_params.clonal_tree_simulator_params.pool_manager_strategy;
    if (pool_manager_strategy == PoolManagerStrategy::UniformPoolManager) {
        return __GetForestStorage<UniformPoolManager>(base_repertoire);
    } else if (pool_manager_strategy == PoolManagerStrategy::DeepTreePoolManager) {
        return __GetForestStorage<DeepTreePoolManager>(base_repertoire);
    } else if (pool_manager_strategy == PoolManagerStrategy::WideTreePoolManager) {
        return __GetForestStorage<WideTreePoolManager>(base_repertoire);
    }
    VERIFY(false);
}

void IgSimulatorLaunch::Run() {
    // MTSingleton::SetSeed(1);
    INFO("== IgSimulator starts ==");

    germline_utils::ChainType chain_type = GetLaunchChainType();
    std::vector<germline_utils::CustomGeneDatabase> db { GetDB(chain_type) };

    const BaseRepertoire base_repertoire = GetBaseRepertoire(chain_type, db);
    const ForestStorage forest_storage = GetForestStorage(base_repertoire);

    INFO("== IgSimulator ends ==");
}

} // End namespace ig_simulator
