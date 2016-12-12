#pragma once

#include "include_me.hpp"

struct BasicRepertoireParams {
    size_t base_repertoire_size;
    size_t mutated_repertoire_size;
    size_t final_repertoire_size;

    BasicRepertoireParams() :
        base_repertoire_size(0),
        mutated_repertoire_size(0),
        final_repertoire_size(0) { }

    BasicRepertoireParams(size_t base_repertoire_size,
        size_t mutated_repertoire_size,
        size_t final_repertoire_size) :
        base_repertoire_size(base_repertoire_size),
        mutated_repertoire_size(mutated_repertoire_size),
        final_repertoire_size(final_repertoire_size) { }
};

ostream& operator<<(ostream &out, const BasicRepertoireParams &params) {
    out << "Base repertoire size: " << params.base_repertoire_size << endl;
    out << "Expected size of mutated repertoire: " << params.mutated_repertoire_size << endl;
    out << "Expected size of final repertoire: " << params.final_repertoire_size << endl;
    return out;
}

struct PatternSHMParams {
    size_t min_number_pattern_shm;
    size_t max_number_pattern_shm;
    double substitution_propability;

    PatternSHMParams(size_t min_number_pattern_shm,
        size_t max_number_pattern_shm,
        double substitution_propability) :
        min_number_pattern_shm(min_number_pattern_shm),
        max_number_pattern_shm(max_number_pattern_shm),
        substitution_propability(substitution_propability) { }

    PatternSHMParams() :
            min_number_pattern_shm(0),
            max_number_pattern_shm(0),
            substitution_propability(0) { }

    static PatternSHMParams CreateStandardParams() {
        return PatternSHMParams(1, 5, .8);
    }
};

struct CDR_SHMParams {
    size_t min_number_mutations;
    size_t max_number_mutations;
    double mutation_in_fr_prop;

    CDR_SHMParams() :
        min_number_mutations(0),
        max_number_mutations(0),
        mutation_in_fr_prop(0) { }

    CDR_SHMParams(size_t min_number_mutations,
                  size_t max_number_mutations,
                  double mutation_in_fr_prop) :
            min_number_mutations(min_number_mutations),
            max_number_mutations(max_number_mutations),
            mutation_in_fr_prop(mutation_in_fr_prop) { }

    static CDR_SHMParams CreateStandardParams() {
        return CDR_SHMParams(1, 5, .25);
    }
};

struct OutputParams {
    string base_sequence_fname;
    string base_multiplicity_fname;
    string mutated_sequence_fname;
    string mutated_multiplicity_fname;
    string mutated_positions;
    string final_repertoire_fname;
    string vdj_recombination_fname;

    OutputParams() :
            base_sequence_fname(),
            base_multiplicity_fname(),
            mutated_sequence_fname(),
            mutated_multiplicity_fname(),
            mutated_positions(),
            final_repertoire_fname(),
            vdj_recombination_fname() { }

    OutputParams(string base_sequence_fname,
        string base_multiplicity_fname,
        string mutated_sequence_fname,
        string mutated_multiplicity_fname,
        string mutated_positions,
        string final_repertoire_fname,
        string vdj_recombination_fname) :
            base_sequence_fname(base_sequence_fname),
            base_multiplicity_fname(base_multiplicity_fname),
            mutated_sequence_fname(mutated_sequence_fname),
            mutated_multiplicity_fname(mutated_multiplicity_fname),
            mutated_positions(mutated_positions),
            final_repertoire_fname(final_repertoire_fname),
            vdj_recombination_fname(vdj_recombination_fname) { }

    static OutputParams CreateStandardParams() {
        return OutputParams("base_sequences.fasta",
                            "base_frequencies.txt",
                            "mutated_sequences.fasta",
                            "mutated_frequencies.txt",
                            "shm_positions.txt",
                            "final_repertoire.fasta",
                            "repertoire_vdj_recombination.txt");
    }

    void AddPrefix(string prefix) {
        base_sequence_fname = prefix + base_sequence_fname;
        base_multiplicity_fname = prefix + base_multiplicity_fname;
        mutated_sequence_fname = prefix + mutated_sequence_fname;
        mutated_multiplicity_fname = prefix + mutated_multiplicity_fname;
        mutated_positions = prefix + mutated_positions;
        final_repertoire_fname = prefix + final_repertoire_fname;
        vdj_recombination_fname = prefix + vdj_recombination_fname;
    }
};

enum DatabaseType {regular_db, imgt_db};

struct HC_InputParams {
    string vgenes_fname;
    string dgenes_fname;
    string jgenes_fname;
    DatabaseType database_type;

    string output_dir;

    OutputParams output_params;
    BasicRepertoireParams basic_repertoire_params;
    PatternSHMParams pattern_shm_params;
    CDR_SHMParams cdr_shm_params;

    HC_InputParams() :
            vgenes_fname(""),
            dgenes_fname(""),
            jgenes_fname(""),
            database_type(regular_db),
            basic_repertoire_params(),
            pattern_shm_params() { }

    void PrintDatabaseParams() {
        cout << "V gene file: " << vgenes_fname << endl;
        cout << "D gene file: " << dgenes_fname << endl;
        cout << "J gene file: " << jgenes_fname << endl;
    }
};

struct LC_InputParams {
    string vgenes_fname;
    string jgenes_fname;
    DatabaseType database_type;

    string output_dir;

    OutputParams output_params;
    BasicRepertoireParams basic_repertoire_params;
    PatternSHMParams pattern_shm_params;
    CDR_SHMParams cdr_shm_params;

    LC_InputParams() :
            vgenes_fname(""),
            jgenes_fname(""),
            database_type(regular_db),
            basic_repertoire_params(),
            pattern_shm_params() { }

    void PrintDatabaseParams() {
        cout << "V gene file: " << vgenes_fname << endl;
        cout << "J gene file: " << jgenes_fname << endl;
    }
};

