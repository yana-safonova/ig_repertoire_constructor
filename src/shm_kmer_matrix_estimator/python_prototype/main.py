#!/usr/bin/env python2

from __future__ import print_function

import argparse
import logging
import os
import sys

from config import Config
from shm_substitution_model_estimator import SHMSubstitutionModelEstimator


def parse_command_line():
    def check_input(file_name):
        if not os.path.isfile(file_name):
            raise argparse.ArgumentTypeError(file_name +
                                             ' does not exists')
        else:
            return file_name

    def prepare_output_directory(output_file):
        directory, filename = os.path.split(output_file)
        if not os.path.exists(directory):
            os.makedirs(directory)
        return output_file

    def check_kmer_length(length):
        if isinstance(length, int) and length > 0 and length & 1:
            return length
        else:
            raise argparse.ArgumentTypeError(length + ' is not a valid \
k-mer length')

    def check_mismatch_finder(mf):
        if mf == 'trivial' or mf == 'k_neighbour':
            return mf
        else:
            raise argparse.ArgumentTypeError(mf + ' is not a valid \
mismatch finder algorithm')

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_filename', type=check_input,
                        help='Alignment input file', required=True)
    parser.add_argument('-o', '--output_filename', type=prepare_output_directory,
                        help='Output file for substitution SHM model',
                        required=True)
    parser.add_argument('-k', '--k_mer_len', type=check_kmer_length,
                        help='K-mer length to estimate the SHM model',
                        default=5)
    parser.add_argument('-mf', '--mismatch_finder', type=check_mismatch_finder,
                        help='What mismatches to take into account: \
                              \'trivial\' / \'k_neighbour\'',
                        required=True)

    command_line_args = parser.parse_args()
    return command_line_args


def prepare_log(config):
    log = logging.getLogger('shm_k-mer_substitution_model')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    log_filename = os.path.join(os.path.dirname(config.IO.output_filename),
                                'shm_k-mer_substitution_model.log')
    log.info(log_filename)
    if os.path.exists(log_filename):
        os.remove(log_filename)
    log_handler = logging.FileHandler(log_filename, mode='a')
    log.addHandler(log_handler)
    log.info('Log will be written to ' + log_filename + '\n')
    return log


def main():
    config = Config(parse_command_line())
    log = prepare_log(config)
    shm_substitution_model_estimator = SHMSubstitutionModelEstimator(config, log)
    shm_substitution_model_estimator.log.info(shm_substitution_model_estimator.config)
    shm_substitution_model_estimator.log.info('========== READING ALIGNMENT STARTED =========')
    shm_substitution_model_estimator.read_alignments()
    shm_substitution_model_estimator.log.info('========== READING ALIGNMENT FINISHED ========')

    shm_substitution_model_estimator.log.info('========== ESTIMATING SUBSTITUTIONS STARTED =========')
    shm_substitution_model_estimator.estimate_substitutions()
    shm_substitution_model_estimator.log.info('========== ESTIMATING SUBSTITUTIONS FINISHED =========')
    shm_substitution_model_estimator.log.info('========== EXPORTING STARTED =========')

    shm_substitution_model_estimator.export_substitutions()
    shm_substitution_model_estimator.log.info('========== EXPORTING FINISHED =========')


if __name__ == '__main__':
    main()
