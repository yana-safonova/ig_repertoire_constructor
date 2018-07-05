#!/usr/bin/env python2

from __future__ import print_function

import numpy as np

import likelihood_calculator.likelihood_calculator as likelihood_calculator
import sample_reader.v_alignments_reader_utilities.alignment_reader as \
    alignment_reader
import shm_kmer_model.shm_kmer_model as shm_kmer_model


class GermlineTester(object):
    def get_likelihood_statistics(self,
                                  model,
                                  model_mode=shm_kmer_model.ModelMode.Both,
                                  input_v_germline_filename='germline' +
                                  '_test_utilities/v_alignment_lymphoma.fasta',
                                  input_v_germline=None,
                                  mismatch_strategy='NoKNeighbours'):
        if input_v_germline is None and input_v_germline_filename is None:
            raise ValueError('No input')
        if input_v_germline is None:
            reader = alignment_reader.GermlineAlignmentReader()
            alignments = reader.read(input_v_germline_filename)
        if input_v_germline_filename is None:
            alignments = input_v_germline

        lkhd_calc = \
            likelihood_calculator.LikelihoodCalculator(model=model,
                                                       model_mode=model_mode)
        results = []
        for ind, alignment in enumerate(alignments):
            print('%d / %d' % (ind, len(alignments)))
            source, dest = alignment.germline_seq.seq, alignment.read.seq
            source, dest = str(source), str(dest)
            results.append((lkhd_calc.calculate_likelihood(source, dest),
                            lkhd_calc.calculate_likelihood(dest, source)))
        return np.array(results, dtype=np.dtype('float, float'))
