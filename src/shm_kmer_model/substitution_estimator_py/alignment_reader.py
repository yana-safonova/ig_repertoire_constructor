#!/usr/bin/env python2

from __future__ import print_function

from Bio import SeqIO

from germline_alignment import GermlineAlignment


class GermlineAlignmentReader:
    def read(self):
        alignments = []
        while True:
            try:
                read = str(next(self.alignment_iterator).seq)
                germline_seq = str(next(self.alignment_iterator).seq)
                alignment = GermlineAlignment(bio_python_read=read,
                                              bio_python_germline=germline_seq)
                if self.alignment_checker.check(alignment):
                    alignments.append(self.alignment_cropper.crop(alignment))
            except StopIteration:
                break
        return alignments

    def __init__(self, input_file, alignment_checker, alignment_checker_params,
                 alignment_cropper, alignment_cropper_params,
                 alignment_format='fasta'):
        self.alignment_iterator = SeqIO.parse(input_file, alignment_format)
        self.alignment_checker = alignment_checker(alignment_checker_params)
        self.alignment_cropper = alignment_cropper(alignment_cropper_params)
