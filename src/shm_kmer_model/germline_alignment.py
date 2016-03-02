#!/usr/bin/env python2

from __future__ import print_function

from Bio import SeqIO


class GermlineAlignment:
    def __init__(self, bio_python_read, bio_python_germline):
        self.read = bio_python_read
        self.germline_seq = bio_python_germline

    def __str__(self):
        return self.read.__str__() + "\n" + self.germline_seq.__str__()

    def __repr__(self):
        return self.read.__repr__() + "\n" + self.germline_seq.__repr__()


class AbstractAlignmentChecker:
    pass


class NoGapsAlignmentChecker(AbstractAlignmentChecker):
    def check(self, alignment, gap_symbol='-'):
        return alignment.read.seq.count(gap_symbol) == 0


class GermlineAlignmentReader:
    def read(self):
        alignments = []
        while True:
            try:
                read = next(self.alignment_iterator)
                germline_seq = next(self.alignment_iterator)
                alignment = GermlineAlignment(bio_python_read=read,
                                              bio_python_germline=germline_seq)
                if self.alignment_checker.check(alignment, self.gap_symbol):
                    alignments.append(alignment)
            except StopIteration:
                break
        return alignments

    def __init__(self, input_file, alignment_checker,
                 alignment_format='fasta', gap_symbol='-',
                ):
        self.alignment_iterator = SeqIO.parse(input_file, alignment_format)
        self.gap_symbol = gap_symbol
        self.alignment_checker = alignment_checker()
