#!/usr/bin/python -u
import logging

from simulation_aux import *
from plot_umi import plot_sens_prec_umi

stages = ("amplify", "barigrec", "igrec", "presto", "migec")
stage_methods = (GetSimulateSteps, GetBarigrecSteps, GetIgrecSteps, GetPrestoSteps, GetMigecSteps)
assert len(stages) == len(stage_methods)
stage_to_method = dict(zip(stages, stage_methods))


def CreateLogger(output_dir):
    log = logging.getLogger('run_sim')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    log_filename = os.path.join(output_dir, "sim.log")
    log_handler = logging.FileHandler(log_filename, mode='a' if os.path.exists(log_filename) else 'w')
    log.addHandler(log_handler)
    log.info("Log will be written to " + log_filename + "\n")
    return log


def ParseCommandLineParams():
    from argparse import ArgumentParser
    parser = ArgumentParser()

    parser.add_argument("-r", "--repertoire",
                        dest = "input_repertoire",
                        type = str,
                        default = "final_repertoire.fasta",
                        help = "Path to input simulated repertoire")
    parser.add_argument("-o", "--output",
                        dest = "output_dir",
                        type = str,
                        default = os.getcwd(),
                        help = "Path to output directory")
    parser.add_argument("-t", "--threads",
                        dest = "threads",
                        type = int,
                        default = 16,
                        help = "Number of threads to use")
    parser.add_argument("-p", "--error-rates",
                        dest = "error_rates_str",
                        type = str,
                        required = True,
                        help = "Comma separated nucleotide error rates")
    parser.add_argument("-b", "--barcode-length",
                        dest = "barcode_length",
                        type = int,
                        default = 15,
                        help = "Length of simulated barcodes")
    parser.add_argument("-e", "--exit_on_error",
                        dest = "exit_on_error",
                        action = "store_true",
                        help = "If set, stop running on first error. By default ignores all errors.")
    stages_args = parser.add_argument_group("Stages arguments. If some of these arguments are specified, only those stages will be run.")
    for stage in stages:
        stages_args.add_argument("--%s" % stage,
                                 dest = stage,
                                 action = "store_true",
                                 help = "Run %s stage" % stage)

    params = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_usage()
        exit(1)

    params.pcr_error_rates = [float(er) for er in params.error_rates_str.split(",")]
    specified_stages = [stage for stage in stages if getattr(params, stage)]
    params.stages_to_run = specified_stages if len(specified_stages) > 0 else stages

    print "Running with params:"
    for param_name in ("input_repertoire", "output_dir", "threads", "pcr_error_rates", "stages_to_run", "exit_on_error"):
        print "%s:%s" % (param_name, getattr(params, param_name))
    print "Working directory: %s" % current_dir

    return params


def RunSimPipeline(run_params, params, log):
    steps = []
    for stage in params.stages_to_run:
        steps.extend(stage_to_method[stage](params, run_params))

    for step in steps:
        exit_status = step.Run(log)
        if exit_status != 0 and params.exit_on_error:
            exit(exit_status)


def main():
    print "Starting run_sum ", sys.argv
    params = ParseCommandLineParams()
    if not os.path.exists(params.output_dir):
        os.makedirs(params.output_dir)
    log = CreateLogger(params.output_dir)
    log.info("Running following stages: %s" % ", ".join(params.stages_to_run))
    log.info("Using error rates: %s" % ", ".join(map(str, params.pcr_error_rates)))

    # pcr_error_rates = [0.006, 0.0006, 0.0025, 0.0012, 0.0018, 0.0030, 0.0036]
    # for supernode_threshold in [100000, 10, 5]:
    for supernode_threshold in [100000]:
        # for barcode_length in [15, 9, 12]:
        for barcode_length in [params.barcode_length]:
            # for pcr_error_rate in [0.006, 0.0006, 0.0025]:
            for pcr_error_rate in params.pcr_error_rates:
                data_path = "%s/pcr_%g_super_%d_umi_%d" % (params.output_dir, pcr_error_rate, supernode_threshold, barcode_length)
                if not os.path.exists(data_path):
                    os.makedirs(data_path)
                simulated_repertoire = os.path.join(data_path, "final_repertoire.fasta")
                if os.path.exists(simulated_repertoire):
                    os.remove(simulated_repertoire)
                os.link(params.input_repertoire, simulated_repertoire)
                RunSimPipeline(RunParams(data_path, pcr_error_rate, supernode_threshold, barcode_length), params, log)

    PyStep("Drawing plots",
           lambda: plot_sens_prec_umi(params.output_dir, params.pcr_error_rates)).Run(log)


if __name__ == '__main__':
    main()
