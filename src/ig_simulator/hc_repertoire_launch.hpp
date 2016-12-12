#include "gene_database.hpp"
#include "repertoire.hpp"
#include "multiplicity_creator.hpp"
#include "config.hpp"
#include "base_repertoire_creation/vdj_recombinator.hpp"
#include "base_repertoire_creation/exonuclease_remover.hpp"
#include "base_repertoire_creation/p_nucleotides_creator.hpp"
#include "base_repertoire_creation/n_nucleotides_creator.hpp"
#include "base_repertoire_creation/cdr_labeler.hpp"
#include "mutated_repertoire_creation/shm_creator.hpp"

ShortGeneNameExtractorPtr CreateShortReadNameExtractor(HC_InputParams params) {
    if(params.database_type == imgt_db)
        return ShortGeneNameExtractorPtr(new IMGTShortGeneNameExtractor());
    return ShortGeneNameExtractorPtr(new TrivialShortGeneNameExtractor());
}

HC_GenesDatabase_Ptr CreateHCDatabase(HC_InputParams params) {
    ShortGeneNameExtractorPtr short_read_extractor_ptr = CreateShortReadNameExtractor(params);
    HC_GenesDatabase_Ptr hc_database(new HC_GenesDatabase(short_read_extractor_ptr));
    hc_database->AddGenesFromFile(variable_gene, params.vgenes_fname);
    hc_database->AddGenesFromFile(diversity_gene, params.dgenes_fname);
    hc_database->AddGenesFromFile(join_gene, params.jgenes_fname);
    return hc_database;
}

HC_Repertoire_Ptr CreateBaseHCRepertoire(HC_InputParams params, HC_GenesDatabase_Ptr hc_database) {
    HC_SimpleRecombinator hc_vdj_recombinator(hc_database, HC_SimpleRecombinationCreator(hc_database));

    // endonuclease remover
    HC_SimpleRemovingStrategy removing_strategy;
    HC_SimpleExonucleaseRemover remover(removing_strategy);

    // p nucleotides insertion
    HC_SimplePInsertionStrategy p_nucls_strategy;
    HC_SimplePNucleotidesCreator p_creator(p_nucls_strategy);

    // n nucleotides insertion
    HC_SimpleNInsertionStrategy n_nucls_strategy;
    HC_SimpleNNucleotidesCreator n_creator(n_nucls_strategy);

    // repertoire
    HC_Repertoire_Ptr base_repertoire(new HC_Repertoire());

    // base multiplicity creator
    //HC_ExponentialMultiplicityCreator base_multiplicity_creator(params.basic_repertoire_params.base_repertoire_size,
    //                                                            params.basic_repertoire_params.mutated_repertoire_size);
    HC_PowerLawMultiplicityCreator base_multiplicity_creator(params.basic_repertoire_params.base_repertoire_size,
                                                          params.basic_repertoire_params.mutated_repertoire_size);

    // cdr labeling
    HC_CDRLabelingStrategy cdr_labeling_strategy;
    HC_CDRLabeler cdr_labeler(cdr_labeling_strategy);

    for(size_t i = 0; i < params.basic_repertoire_params.base_repertoire_size; i++) {
        auto vdj = hc_vdj_recombinator.CreateRecombination();
        vdj = remover.CreateRemovingSettings(vdj);
        vdj = p_creator.CreatePNucleotides(vdj);
        vdj = n_creator.CreateNNucleotides(vdj);

        HC_VariableRegionPtr ig_variable_region = HC_VariableRegionPtr(new HC_VariableRegion(vdj));
        ig_variable_region = cdr_labeler.LabelCDRs(ig_variable_region);

        size_t multiplicity = base_multiplicity_creator.AssignMultiplicity(ig_variable_region);
        base_repertoire->Add(HC_Cluster(ig_variable_region, multiplicity));
    }
    return base_repertoire;
}

void PrintBaseHCRepertoire(HC_Repertoire_Ptr base_repertoire) {
    for(auto it = base_repertoire->begin(); it != base_repertoire->end(); it++) {
        cout << it->IgVariableRegion()->VDJ_Recombination()->Sequence() << " " << it->Multiplicity() << endl;
        cout << it->IgVariableRegion()->GetCDRSettings();
    }
}

HC_Repertoire_Ptr CreateMutatedHCRepertoire(HC_InputParams params, HC_Repertoire_Ptr base_repertoire) {
    HC_Repertoire_Ptr mutated_repertoire(new HC_Repertoire());

    // mutated multiplicity creator
//    HC_ExponentialMultiplicityCreator mutated_multiplicity_creator(base_repertoire->NumberAntibodies(),
//        params.basic_repertoire_params.final_repertoire_size);
    HC_PowerLawMultiplicityCreator mutated_multiplicity_creator(base_repertoire->NumberAntibodies(),
                                                          params.basic_repertoire_params.final_repertoire_size);

    // shm creator
    HC_RgywWrcySHMStrategy shm_creation_strategy1(params.pattern_shm_params.min_number_pattern_shm,
                                                 params.pattern_shm_params.max_number_pattern_shm,
                                                 params.pattern_shm_params.substitution_propability);
    HC_CDRBasedRandomSHMStrategy shm_creation_strategy2(params.cdr_shm_params.min_number_mutations,
                                                        params.cdr_shm_params.max_number_mutations,
                                                        params.cdr_shm_params.mutation_in_fr_prop);
    HC_CompositeSHMCreationStrategy composite_shm_strategy(shm_creation_strategy1, shm_creation_strategy2);
    HC_SHMCreator shm_creator(composite_shm_strategy);

    for(auto it = base_repertoire->begin(); it != base_repertoire->end(); it++) {
        for(size_t i = 0; i < it->Multiplicity(); i++) {
            auto variable_region_ptr = it->IgVariableRegion()->Clone();
            variable_region_ptr = shm_creator.CreateSHM(variable_region_ptr);
            size_t multiplicity = mutated_multiplicity_creator.AssignMultiplicity(variable_region_ptr);
            mutated_repertoire->Add(HC_Cluster(variable_region_ptr, multiplicity));
        }
    }
    return mutated_repertoire;
}

void PrintMutatedHCRepertoire(HC_Repertoire_Ptr mutated_repertoire) {
    for(auto it = mutated_repertoire->begin(); it != mutated_repertoire->end(); it++) {
        cout << it->IgVariableRegion()->VDJ_Recombination()->Sequence() << " " << it->Multiplicity() << endl;
    }
}

void CreateHCRepertoire(HC_InputParams params) {
    cout << "======== Simulation of heavy chain repertoire starts" << endl;
    cout << "Repertoire parameters: " << endl;
    cout << params.basic_repertoire_params << endl;

    cout << "Database parameters: " << endl;
    params.PrintDatabaseParams();
    cout << endl;

    auto hc_database = CreateHCDatabase(params);

    cout << "==== Generation of base repertoire" << endl;
    HC_Repertoire_Ptr base_repertoire = CreateBaseHCRepertoire(params, hc_database);
    cout << "Base repertoire consists of " << base_repertoire->Size() << " sequences with total multiplicities " <<
            base_repertoire->NumberAntibodies() << endl;
    base_repertoire->OutputSequences(params.output_params.base_sequence_fname);
    base_repertoire->OutputMultiplicities(params.output_params.base_multiplicity_fname);
    cout << "Base antibody sequences were written to " << params.output_params.base_sequence_fname << endl;
    cout << "Base antibody multiplicities were written to " << params.output_params.base_multiplicity_fname << endl;
    cout << endl;

    cout << "==== Generation of mutated repertoire" << endl;
    HC_Repertoire_Ptr mutated_repertoire = CreateMutatedHCRepertoire(params, base_repertoire);
    cout << "Mutated repertoire consists of " << mutated_repertoire->Size() <<
            " sequences with total multiplicities " << mutated_repertoire->NumberAntibodies() << endl;
    mutated_repertoire->OutputSequences(params.output_params.mutated_sequence_fname);
    mutated_repertoire->OutputMultiplicities(params.output_params.mutated_multiplicity_fname);
    mutated_repertoire->OutputSHMPositions(params.output_params.mutated_positions);
    cout << "Mutated antibody sequences were written to " << params.output_params.mutated_sequence_fname << endl;
    cout << "Mutated antibody multiplicities were written to " <<
            params.output_params.mutated_multiplicity_fname << endl;
    cout << "Positions of SHM were written to " << params.output_params.mutated_positions << endl;
    mutated_repertoire->OutputVDJRecombination(params.output_params.vdj_recombination_fname);
    cout << "VDJ recombination of the final repertoire was written to " <<
    params.output_params.vdj_recombination_fname << endl;
    cout << endl;

    mutated_repertoire->OutputRepertoire(params.output_params.final_repertoire_fname);
    cout << "Final repertoire was written to " << params.output_params.final_repertoire_fname << endl << endl;

    cout << "======== Simulation of heavy chain repertoire ends" << endl;
}