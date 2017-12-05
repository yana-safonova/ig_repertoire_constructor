#!/usr/bin/env python2

import os
import sys
import init
import logging
import shutil
import ntpath

home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + '/'
py_src = os.path.join(home_directory, "src/python_pipeline/")

antevolo_config_dir = os.path.join(home_directory, "configs/antevolo/")
cdr_labeler_config_dir = os.path.join(home_directory, "configs/cdr_labeler")
vj_finder_config_dir = os.path.join(home_directory, "configs/vj_finder")
shm_kmer_matrix_estimator_config_dir = os.path.join(home_directory, "configs/shm_kmer_matrix_estimator")

antevolo_bin = os.path.join(home_directory, "build/release/bin/antevolo")
run_antevolo = os.path.join(home_directory, "build/release/bin/./antevolo")
hg_constructor_bin = os.path.join(home_directory, "build/release/bin/ig_swgraph_construct")
run_hg_constructor = os.path.join(home_directory, "build/release/bin/./ig_swgraph_construct")

sys.path.append(py_src)
#sys.path.append(visualizer_dir)
import process_cfg
import support

test_reads = os.path.join(home_directory, "test_dataset/merged_reads.fastq")
test_dir = os.path.join(os.getcwd(), "antevolo_test")

tool_name = "AntEvolo"
comparing_mode_name = "EvoQuast"

def CheckBinariesExistance(params, log):
    if not os.path.exists(antevolo_bin):
        log.info("ERROR: " + tool_name + " binary file was not found. Please type \'make\' before running.")
        sys.exit(1)
    if not os.path.exists(hg_constructor_bin):
        log.info("ERROR: HG constructor binary file was not found. Please type \'make\' before running.")

def LociParamCorrect(loci_str):
    return loci_str == "IG" or loci_str == "IGH" or loci_str == "IGK" or loci_str == "IGL"

def CheckParamsCorrectness(params, log):
    if not os.path.exists(params.input_reads):
        log.info("Input reads " + params.input_reads + " were not found")
        sys.exit(1)
    if not LociParamCorrect(params.loci):
        log.info("Loci " + params.loci + " is not recognized")
        sys.exit(1)
    if params.compare and not os.path.exists(params.rcm):
        log.info("No RCM in " + comparing_mode_name + " mode: RCM file " + params.rcm + " was not found")
        sys.exit(1)

def SetOutputParams(params, log):
    if params.input_reads == test_reads:
        params.output_dir = test_dir
    if params.input_reads != test_reads and params.output_dir == "":
        log.info("ERROR: Output dir (-o) was not specified")
        sys.exit(1)
    params.output_dir = os.path.abspath(params.output_dir)
    params.config_dir = os.path.join(params.output_dir, "configs")
    params.antevolo_config_file = os.path.join(antevolo_config_dir, "config.info")
    params.cdr_labeler_config_file = os.path.join(cdr_labeler_config_dir, "config.info")
    params.vj_finder_config_file = os.path.join(vj_finder_config_dir, "config.info")

def PrepareOutputDir(params):
    if os.path.exists(params.output_dir):
        shutil.rmtree(params.output_dir)
    os.mkdir(params.output_dir)

def PrintParams(params, log):
    log.info(tool_name + " parameters:")
    log.info("  Mode:\t" + (comparing_mode_name if params.compare else "default"))
    log.info("  Input reads:\t\t\t" + params.input_reads)
    if params.compare:
        log.info("  Input decomposition:\t\t" + params.rcm)
    log.info("  Output directory:\t\t" + params.output_dir)
    log.info("  Number of threads:\t\t" + str(params.num_threads) + "\n")
    log.info("  Loci:\t\t\t\t" + params.loci)
    log.info("  Max distance between CDR3:\t" + str(params.cdr3_tau))
    log.info("  Min number of shared SHMs:\t" + str(params.min_shared_shms))
    log.info("  Use SHM model:\t\t" + str(params.model))


########################################################################################################################

def CopyConfigs(params, log):
    if os.path.exists(params.config_dir):
        shutil.rmtree(params.config_dir)
    os.mkdir(params.config_dir)
    # copying antevolo config
    params.antevolo_config_dir = os.path.abspath(os.path.join(params.config_dir, "antevolo"))
    shutil.copytree(antevolo_config_dir, params.antevolo_config_dir)
    # copying cdr labeler config
    params.cdr_labeler_config_dir = os.path.abspath(os.path.join(params.config_dir, "cdr_labeler"))
    shutil.copytree(cdr_labeler_config_dir, params.cdr_labeler_config_dir)
    # copying vj finder config
    params.vj_finder_config_dir = os.path.abspath(os.path.join(params.config_dir, "vj_finder"))
    shutil.copytree(vj_finder_config_dir, params.vj_finder_config_dir)
    # copying shm kmer matrix estimator config
    params.shm_kmer_matrix_estimator_config_dir = os.path.abspath(os.path.join(params.config_dir, "shm_kmer_matrix_estimator"))
    shutil.copytree(shm_kmer_matrix_estimator_config_dir, params.shm_kmer_matrix_estimator_config_dir)
    # checking config files
    params.antevolo_config_file = os.path.join(params.antevolo_config_dir, "config.info")
    params.cdr_labeler_config_file = os.path.join(params.cdr_labeler_config_dir, "config.info")
    params.vj_finder_config_file = os.path.join(params.vj_finder_config_dir, "config.info")
    params.shm_kmer_matrix_estimator_config_file = os.path.join(params.shm_kmer_matrix_estimator_config_dir, "config.info") # s?? wtf?
    params.germline_config_file = os.path.join(params.vj_finder_config_dir, "germline_files_config.txt")

    if not os.path.exists(params.antevolo_config_dir):
        log.info("ERROR: Config file " + params.antevolo_config_dir + " was not found")
        sys.exit(1)
    if not os.path.exists(params.vj_finder_config_file):
        log.info("ERROR: Config file " + params.vj_finder_config_file + " was not found")
        sys.exit(1)
    if not os.path.exists(params.cdr_labeler_config_file):
        log.info("ERROR: Config file " + params.cdr_labeler_config_file + " was not found")
        sys.exit(1)
    if not os.path.exists(params.shm_kmer_matrix_estimator_config_file):
        log.info("ERROR: Config file " + params.shm_kmer_matrix_estimator_config_file + " was not found")
        sys.exit(1)

def UpdateGermlineConfigFile(params, log):
    lines = open(params.germline_config_file, "r").readlines()
    fhandler = open(params.germline_config_file, "w")
    fhandler.write(lines[0].strip() + "\n")
    for i in range(1, len(lines)):
        splits = lines[i].strip().split()
        fhandler.write(splits[0] + "\t" + splits[1] + "\t" + os.path.join(home_directory, splits[2]) + "\n")
    fhandler.close()

def ModifyParallelModeParams(params, param_dict):
    if not params.parallel_survey:
        return
    param_dict['min_num_intersected_v_shms'] = 400
    param_dict['enable_parallel_shms_finder'] = 'true'

def ModifyAntEvoloConfigFile(params, log):
    param_dict = dict()
    param_dict['input_reads'] = params.input_reads
    param_dict['decomposition_rcm'] = params.rcm
    param_dict['output_dir'] = params.output_dir
    param_dict['cdr_labeler_config_fname'] = params.cdr_labeler_config_file
    param_dict['compare'] = int(params.compare)
    param_dict['model'] = int(params.model)
    param_dict['shm_kmer_matrix_estimator_config_fname'] = params.shm_kmer_matrix_estimator_config_file
    param_dict['shm_kmer_model_igh'] = os.path.join(home_directory, "data/shm_model/NoKNeighbours_IGH.csv")
    param_dict['shm_kmer_model_igk'] = os.path.join(home_directory, "data/shm_model/NoKNeighbours_IGK.csv")
    param_dict['shm_kmer_model_igl'] = os.path.join(home_directory, "data/shm_model/NoKNeighbours_IGL.csv")
    param_dict['num_threads'] = params.num_threads
    ModifyParallelModeParams(params, param_dict)
    process_cfg.substitute_params(params.antevolo_config_file, param_dict, log)

def ModifyCDRLabelerConfigFile(params, log):
    param_dict = dict()
    param_dict['num_threads'] = params.num_threads
    param_dict['output_dir'] = params.output_dir
    param_dict['vj_finder_config'] = params.vj_finder_config_file
    param_dict['imgt_v_annotation'] = os.path.join(home_directory, 'data/annotation/human_v_imgt.txt')
    param_dict['kabat_v_annotation'] = os.path.join(home_directory, 'data/annotation/human_v_kabat.txt')
    param_dict['imgt_j_annotation'] = os.path.join(home_directory, 'data/annotation/human_j_imgt.txt')
    param_dict['kabat_j_annotation'] = os.path.join(home_directory, 'data/annotation/human_j_kabat.txt')
    param_dict['run_hg_constructor'] = os.path.join(home_directory, './build/release/bin/ig_swgraph_construct')
    process_cfg.substitute_params(params.cdr_labeler_config_file, param_dict, log)

def ModifyVjFinderConfigFile(params, log):
    vj_param_dict = dict()
    vj_param_dict['loci'] = params.loci
    vj_param_dict['germline_dir'] = os.path.join(home_directory, "data/antevolo_germline")
    params.germline_config_file = os.path.join(params.vj_finder_config_dir, "germline_files_config.txt")
    vj_param_dict['germline_filenames_config'] = params.germline_config_file
    process_cfg.substitute_params(params.vj_finder_config_file, vj_param_dict, log)

def ModifyConfigFiles(params, log):
    ModifyAntEvoloConfigFile(params, log)
    ModifyCDRLabelerConfigFile(params, log)
    ModifyVjFinderConfigFile(params, log)

def PrepareConfigs(params, log):
    CopyConfigs(params, log)
    ModifyConfigFiles(params, log)
    #UpdateGermlineConfigFile(params, log)

########################################################################################################################

def main(argv):
    from argparse import ArgumentParser
    parser = ArgumentParser(description="== " + tool_name + ": clonal tree construction algorithm ==",
                                  epilog="In case you have troubles running " + tool_name +
                                         ", you can write to igtools_support@googlegroups.com. " +
                                         "Please provide us with antevolo.log file from the output directory.",
                                  add_help=False)
    req_args = parser.add_argument_group("Required params")
    input_args = req_args.add_mutually_exclusive_group(required=True)
    input_args.add_argument("-i", "--input",
                            type=str,
                            default="",
                            dest="input_reads",
                            help="Input reads in FASTQ/FATSA format")
    input_args.add_argument("--test",
                            action="store_const",
                            const=test_reads,
                            dest="input_reads",
                            help="Running test dataset")

    out_args = parser.add_argument_group("Output")
    out_args.add_argument("-o", "--output",
                          type=str,
                          dest="output_dir",
                          default="", #os.path.join(home_directory, "cdr_test"),
                          help="Output directory")

    optional_args = parser.add_argument_group("Optional arguments")
    optional_args.add_argument("-t", "--threads",
                               type=int,
                               default=16,
                               dest="num_threads",
                               help="Threads number [default: %(default)d]")
    optional_args.add_argument("-l", "--loci",
                               type=str,
                               default="IG",
                               dest="loci",
                               help="Loci: IGH, IGK, IGL, IG (all BCRs)" # ", TRA, TRB, TRG, TRD, TR (all TCRs) or all. "
                                    "[default: %(default)s]")
    optional_args.add_argument("--compare",
                               action="store_true",
                               dest="compare",
                               help="Run " + comparing_mode_name + " to assess the quality of clonal lineages decomposition " +
                                    "instead of reconstructing the trees. Is this case, " +
                                    "RCM file with the decomposition have to be provided")
    optional_args.add_argument("--rcm", "-R",
                               type=str,
                               default="",
                               help="RCM file with the decomposition")
    optional_args.set_defaults(compare=False)

    algorithm_args = parser.add_argument_group("Algorithm arguments")
    algorithm_args.add_argument("--tau",
                                type=int,
                                default=3,
                                dest="cdr3_tau",
                                help="Threshold for construction Hamming graph on CDR3s [default: %(default)d]")
    algorithm_args.add_argument("-s", "--min-shared-shms",
                                type=int,
                                default=3,
                                dest="min_shared_shms",
                                help="Minimal number of SHMs for considering two antibodies clonally related")
    algorithm_args.add_argument("--model",
                                action="store_true",
                                dest="model",
                                help="Use SHM statistical model to improve trees reconstruction")
    algorithm_args.set_defaults(model=False)
    algorithm_args.add_argument("--par-shms",
                                action='store_true',
                                dest='parallel_survey',
                                help='Enables analysis of parallel SHMs, disables reconstruction of missing clones')

    optional_args.add_argument("-h", "--help",
                               action="help",
                               help="Help message and exit")

    #parser.set_defaults(config_dir=os.path.join(params.output_dir, "configs"),
    #                    cdr_config_file=os.path.join(cdr_labeler_config_dir, "config.info"),
    #                    vj_finder_config_file = os.path.join(vj_finder_config_dir, "config.info")))

    # prepare log
    log = logging.getLogger('antevolo')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    params = parser.parse_args()

    CheckBinariesExistance(params, log)
    CheckParamsCorrectness(params, log)
    SetOutputParams(params, log)

    PrepareOutputDir(params)

    # log file
    params.log_filename = os.path.join(params.output_dir, "antevolo.log")
    if os.path.exists(params.log_filename):
        os.remove(params.log_filename)
    log_handler = logging.FileHandler(params.log_filename, mode='a')
    log.addHandler(log_handler)

    # print command line
    command_line = "Command_line: "
    command_line += " ".join(argv)
    log.info(command_line + "\n")
    PrintParams(params, log)
    log.info("Log will be written to " + params.log_filename + "\n")

    PrepareConfigs(params, log)
    try:
        antevolo_command_line = run_antevolo + " " + params.antevolo_config_file
        support.sys_call(antevolo_command_line, log)
        #Cleanup(params, log)
        log.info("\nThank you for using " + tool_name + "!\n")
    except (KeyboardInterrupt):
         log.info("\n" + tool_name + " was interrupted!")
    except Exception:
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            log.exception(exc_value)
            log.info("\nERROR: Exception caught.")
             #supportInfo(log)
    except BaseException:
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            log.exception(exc_value)
            log.info("\nERROR: Exception caught.")
            #supportInfo(log)

    log.info("Log was written to " + params.log_filename)


if __name__ == '__main__':
    main(sys.argv)
