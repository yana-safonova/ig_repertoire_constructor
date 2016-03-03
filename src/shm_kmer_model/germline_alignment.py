#!/usr/bin/env python2

from __future__ import print_function

from Bio import SeqIO


class GermlineAlignment:
    def __init__(self, bio_python_read, bio_python_germline):
        self.read = bio_python_read
        self.germline_seq = bio_python_germline

    def __str__(self):
        result = "Alignment to germline:\n"
        result += "Read:\n"
        result += "\t" + self.read.__str__() + "\n"
        result += "Germline:\n"
        result += "\t" + self.germline_seq.__str__() + "\n"
        return result

    def __repr__(self):
        result = "GermlineAlignment("
        result += "Read="
        result += self.read.__repr__()
        result += ", Germline="
        result += self.germline_seq.__repr__() + ")"
        return result


class AbstractAlignmentChecker:
    pass


class NoGapsNoNAlignmentChecker(AbstractAlignmentChecker):
    # TODO without germline check after @eodus answer
    def check(self, alignment):
        return alignment.read.seq.count(self.gap_symbol) == 0 and \
            alignment.read.seq.count('N') == 0 and \
            alignment.germline_seq.seq.count(self.gap_symbol) == 0 and \
            alignment.germline_seq.seq.count('N') == 0

    def __str__(self):
        return "Alignment Checker: only sequences with no gaps are considered\n"

    def __init__(self, config):
        self.gap_symbol = config.gap_symbol


class AbstractAlignmentCropper:
    pass


class UptoLastReliableKMerAlignmentCropper(AbstractAlignmentCropper):
    def crop(self, alignment):
        """ Crop the tail """
        start_last_reliable_k_mer = alignment.read.__len__() - self.k_mer_len
        while alignment.read.seq[start_last_reliable_k_mer:
                                 start_last_reliable_k_mer + self.k_mer_len] != \
                alignment.germline_seq.seq[start_last_reliable_k_mer:
                                           start_last_reliable_k_mer + self.k_mer_len] and \
                start_last_reliable_k_mer >= 0:
            start_last_reliable_k_mer -= 1

        alignment.read.seq = alignment.read.seq[:start_last_reliable_k_mer + self.k_mer_len]
        alignment.germline_seq.seq = alignment.germline_seq.seq[:start_last_reliable_k_mer +
                                                                self.k_mer_len]

        """ Crop the head """
        start_first_reliable_k_mer = 0
        while alignment.read.seq[start_first_reliable_k_mer:
                                 start_first_reliable_k_mer + self.k_mer_len] != \
                alignment.germline_seq.seq[start_first_reliable_k_mer:
                                           start_first_reliable_k_mer + self.k_mer_len] and \
                start_first_reliable_k_mer >= 0:
            start_first_reliable_k_mer += 1

        alignment.read.seq = alignment.read.seq[start_first_reliable_k_mer:]
        alignment.germline_seq.seq = alignment.germline_seq.seq[start_first_reliable_k_mer:]

        return alignment

    def __init__(self, config):
        self.k_mer_len = config.k_mer_len


class GermlineAlignmentReader:
    def read(self):
        alignments = []
        while True:
            try:
                read = next(self.alignment_iterator)
                germline_seq = next(self.alignment_iterator)
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
