#!/usr/bin/env python2

import os
import sys

current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = current_dir + "/../"
sys.path.append(igrec_dir)
import aimquast


if __name__ == "__main__":
    ig_simulator_output_dir = "/tmp/ig_simulator"
    output_dir = igrec_dir + "/aimquast_test_dataset"
    aimquast.run_ig_simulator(ig_simulator_output_dir)
    aimquast.simulate_data(ig_simulator_output_dir + "/final_repertoire.fasta", output_dir)

    aimquast.run_igrec(output_dir + "/merged_reads.fq",
                       output_dir + "/igrec_good/")

    aimquast.run_igrec(output_dir + "/merged_reads.fq",
                       output_dir + "/igrec_bad/", additional_args="--create-triv-dec")
