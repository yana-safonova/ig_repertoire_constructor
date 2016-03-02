#!/usr/bin/env python2

from __future__ import print_function
import os

from Bio import SeqIO

from germline_alignment import GermlineAlignmentReader


class SHMSubstitutionModelEstimator:
    def read_alignments(self):
       self.alignments = GermlineAlignmentReader(input_file=self.config.IO.input_filename,
                                                 alignment_checker=
                                                 self.config.alignment_checker).read()
       return self.alignments

    def __init__(self, config, log):
        self.config = config
        self.log = log
        self.log.info(config)
