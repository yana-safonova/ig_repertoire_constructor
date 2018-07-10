#!/usr/bin/env python2

from __future__ import print_function

from Bio import SeqIO

from germline_alignment import GermlineAlignment
from alignment_cropper import UptoLastReliableKMerAlignmentCropper
from alignment_checker import NoGapsAlignmentChecker


class GermlineAlignmentReader(object):
    def read(self, input_file):
        alignment_iterator = SeqIO.parse(input_file, self.alignment_format)
        alignments = []
        while True:
            try:
                read = next(alignment_iterator)
                germline_seq = next(alignment_iterator)
                alignment = GermlineAlignment(bio_python_read=read,
                                              bio_python_germline=germline_seq)
                if self.alignment_checker.check(alignment):
                    alignments.append(self.alignment_cropper.crop(alignment))
            except StopIteration:
                break
        return alignments

    def __init__(self,
                 checker=NoGapsAlignmentChecker(),
                 cropper=UptoLastReliableKMerAlignmentCropper(),
                 alignment_format='fasta'):
        self.alignment_checker = checker
        self.alignment_cropper = cropper
        self.alignment_format = alignment_format
