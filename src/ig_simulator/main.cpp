#include "config.hpp"
#include "hc_repertoire_launch.hpp"
#include "lc_repertoire_launch.hpp"


/*
 * ./ig_simulator HC output_dir base_rep_size mutated_rep_size final_rep_size Vgene.fa Dgene.fa Jgene.fa
 * ./ig_simulator LC output_dir base_rep_size mutated_rep_size final_rep_size Vgene.fa Jgene.fa
 */

void HCUsage() {
    cout << "Usage for simulation heavy chain repertoire:" << endl;
    cout << "./ig_simulator HC output_dir base_rep_size mutated_rep_size final_rep_size Vgene.fa Dgene.fa Jgene.fa" << endl;
}

void LCUsage() {
    cout << "Usage for simulation light chain repertoire:" << endl;
    cout << "./ig_simulator LC output_dir base_rep_size mutated_rep_size final_rep_size Vgene.fa Jgene.fa" << endl;
}

void Usage() {
    HCUsage();
    LCUsage();
}

enum ChainType {Unknown_chain, Heavy_chain, Light_chain};

ChainType ParseChainType(string arg) {
    if(arg == "HC")
        return Heavy_chain;
    else if(arg == "LC")
        return Light_chain;
    else
        return Unknown_chain;
}

void PrepareOutputDir(string output_dir) {
    struct stat st = {0};
    if (stat(output_dir.c_str(), &st) == -1) {
        mkdir(output_dir.c_str(), 0700);
    }
}

struct HCParamsIndices {
    static const int output_dir_ind = 2;
    static const int base_size_ind = 3;
    static const int mutated_size_ind = 4;
    static const int final_size_ind = 5;
    static const int vgenes_ind = 6;
    static const int dgenes_ind = 7;
    static const int jgenes_ind = 8;
    static const int database_type = 9;
    static const int num_params = 10;
};

DatabaseType ParseDatabaseType(string database_type_str) {
    if(database_type_str == "imgt")
        return imgt_db;
    return regular_db;
}

HC_InputParams ParseHCInputParams(int argc, char* argv[]) {
    HC_InputParams input_params;
    input_params.output_dir = string(argv[HCParamsIndices::output_dir_ind]);
    if(input_params.output_dir[input_params.output_dir.size() - 1] != '/')
        input_params.output_dir += "/";
    input_params.basic_repertoire_params.base_repertoire_size = StringToType<int>(
            string(argv[HCParamsIndices::base_size_ind]));
    input_params.basic_repertoire_params.mutated_repertoire_size = StringToType<int>(
            string(argv[HCParamsIndices::mutated_size_ind]));
    input_params.basic_repertoire_params.final_repertoire_size = StringToType<int>(
            string(argv[HCParamsIndices::final_size_ind]));
    input_params.vgenes_fname =  string(argv[HCParamsIndices::vgenes_ind]);
    input_params.dgenes_fname =  string(argv[HCParamsIndices::dgenes_ind]);
    input_params.jgenes_fname =  string(argv[HCParamsIndices::jgenes_ind]);
    input_params.database_type = ParseDatabaseType(string(argv[HCParamsIndices::database_type]));

    input_params.output_params = OutputParams::CreateStandardParams();
    input_params.output_params.AddPrefix(input_params.output_dir);

    input_params.pattern_shm_params = PatternSHMParams::CreateStandardParams();
    input_params.cdr_shm_params = CDR_SHMParams::CreateStandardParams();
    return input_params;
}

struct LCParamsIndices {
    static const int output_dir_ind = 2;
    static const int base_size_ind = 3;
    static const int mutated_size_ind = 4;
    static const int final_size_ind = 5;
    static const int vgenes_ind = 6;
    static const int jgenes_ind = 7;
    static const int database_type = 8;
    static const int num_params = 9;
};

LC_InputParams ParseLCInputParams(int argc, char* argv[]) {
    LC_InputParams input_params;
    input_params.output_dir = string(argv[LCParamsIndices::output_dir_ind]);
    if(input_params.output_dir[input_params.output_dir.size() - 1] != '/')
        input_params.output_dir += "/";
    input_params.basic_repertoire_params.base_repertoire_size = StringToType<int>(
            string(argv[LCParamsIndices::base_size_ind]));
    input_params.basic_repertoire_params.mutated_repertoire_size = StringToType<int>(
            string(argv[LCParamsIndices::mutated_size_ind]));
    input_params.basic_repertoire_params.final_repertoire_size = StringToType<int>(
            string(argv[LCParamsIndices::final_size_ind]));
    input_params.vgenes_fname =  string(argv[LCParamsIndices::vgenes_ind]);
    input_params.jgenes_fname =  string(argv[LCParamsIndices::jgenes_ind]);
    input_params.database_type = ParseDatabaseType(string(argv[LCParamsIndices::database_type]));

    input_params.output_params = OutputParams::CreateStandardParams();
    input_params.output_params.AddPrefix(input_params.output_dir);

    input_params.pattern_shm_params = PatternSHMParams::CreateStandardParams();
    input_params.cdr_shm_params = CDR_SHMParams::CreateStandardParams();
    return input_params;
}

int main(int argc, char *argv[]) {
    if(argc < 3) {
        cout << "ERROR: Invalid number of input parameters" << endl;
        Usage();
        return 1;
    }

    string chain_type_arg = string(argv[1]);
    ChainType chain_type = ParseChainType(chain_type_arg);
    if(chain_type == Unknown_chain) {
        cout << "ERROR: Invalid chain type: " << chain_type_arg << endl;
        Usage();
        return 1;
    }

    string output_dir = string(argv[2]);
    PrepareOutputDir(output_dir);
    cout << "Repertoire and statistics will written be to " << output_dir << endl << endl;

    if(chain_type == Heavy_chain) {
        if(argc != HCParamsIndices::num_params) {
            cout << "ERROR: Invalid number of input parameters" << endl;
            HCUsage();
            return 1;
        }
        HC_InputParams input_params = ParseHCInputParams(argc, argv);
        CreateHCRepertoire(input_params);
    }
    else {
        if(argc != LCParamsIndices::num_params) {
            cout << "ERROR: Invalid number of input parameters" << endl;
            LCUsage();
            return 1;
        }
        LC_InputParams input_params = ParseLCInputParams(argc, argv);
        CreateLCRepertoire(input_params);
    }
}