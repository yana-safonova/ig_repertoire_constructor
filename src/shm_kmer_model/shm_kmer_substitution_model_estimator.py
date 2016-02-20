#!/usr/bin/env python2

from __future__ import print_function
import argparse
import logging
import os
import sys


class SHM_KmerSubstitutionModelEstimator:
    def __parse_command_line__(self):
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
                raise argparse.ArgumentTypeError(length + ' is not a valid\
                        k-mer length')

        parser = argparse.ArgumentParser()
        parser.add_argument('-i', '--alignment_input', type=check_input,
                            help='Alignment input file', required=True)
        parser.add_argument('-o', '--output', type=prepare_output_directory,
                            help='Output file for substitution SHM model',
                            required=True)
        parser.add_argument('-k', '--k-mer_length', type=check_kmer_length,
                            help='K-mer length to estimate the SHM model',
                            default=5)
        command_line_args = parser.parse_args()
        return command_line_args

    def __prepare_log__(self, params):
        log = logging.getLogger('shm_k-mer_substitution_model')
        log.setLevel(logging.DEBUG)
        console = logging.StreamHandler(sys.stdout)
        console.setFormatter(logging.Formatter('%(message)s'))
        console.setLevel(logging.DEBUG)
        log.addHandler(console)

        log_filename = os.path.join(os.path.dirname(params.output),
                                    "shm_k-mer_substitution_model.log")
        log.info(log_filename)
        if os.path.exists(log_filename):
            os.remove(log_filename)
        log_handler = logging.FileHandler(log_filename, mode='a')
        log.addHandler(log_handler)
        log.info("Log will be written to " + log_filename + "\n")
        return log

    def __init__(self):
        params = self.__parse_command_line__()
        log = self.__prepare_log__(params)
        log.info(params)
