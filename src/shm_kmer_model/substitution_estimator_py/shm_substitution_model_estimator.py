#!/usr/bin/env python2

from __future__ import print_function
import itertools

import pandas as pd

from alignment_reader import GermlineAlignmentReader


class SHMSubstitutionModelEstimator:
    def read_alignments(self):
        self.alignments = GermlineAlignmentReader(input_file=self.config.IO.input_filename,
                                                  alignment_checker=
                                                  self.config.alignment_checker.method,
                                                  alignment_checker_params=
                                                  self.config.alignment_checker.params,
                                                  alignment_cropper=
                                                  self.config.alignment_cropper.method,
                                                  alignment_cropper_params=
                                                  self.config.alignment_cropper.params).read()
        return self.alignments

    def estimate_substitutions(self):
        bases = ['A', 'C', 'G', 'T']  # TODO Biopython
        all_k_mers = [''.join(p) for p in itertools.product(bases, repeat=self.config.k_mer_len)]
        self.substitution_dataframe = pd.DataFrame(index=all_k_mers, columns=bases)
        self.substitution_dataframe = self.substitution_dataframe.fillna(0)

        mismatch_finder = self.config.mismatch_finder.method(
                config=self.config.mismatch_finder.params)

        modulo_log = 500
        for idx, alignment in enumerate(self.alignments):
            if idx % modulo_log == 0:
                self.log.info('Alignments: processed %d / %d' % (idx, len(self.alignments)))
            mismatch_positions = mismatch_finder.find_mismatch_positions(alignment.read.seq,
                                                                   alignment.germline_seq.seq)
            for mismatch_position in mismatch_positions:
                k_mer = alignment.germline_seq.seq[mismatch_position - self.config.k_mer_len // 2:
                                                   mismatch_position + self.config.k_mer_len // 2 + 1]
                k_mer = str(k_mer)
                if 'N' not in k_mer and 'N' != alignment.read.seq[mismatch_position]:
                    self.substitution_dataframe.ix[k_mer, alignment.read.seq[mismatch_position]] += 1
        return self.substitution_dataframe

    def export_substitutions(self):
        self.substitution_dataframe.to_csv(self.config.IO.output_filename, sep=',')

    def __init__(self, config, log):
        self.config = config
        self.log = log
