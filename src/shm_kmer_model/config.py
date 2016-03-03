#!/usr/bin/env python2

from __future__ import print_function

from germline_alignment import NoGapsNoNAlignmentChecker
from germline_alignment import UptoLastReliableKMerAlignmentCropper


# TODO from config.file
class Config:
    class IOConfig:
        def __init__(self, command_line_args):
            self.input_filename = command_line_args.input_filename
            self.output_filename = command_line_args.output_filename

        def __str__(self):
            return (('Input filename:  %s\n'
                     'Output filename: %s\n') %
                    (self.input_filename,
                     self.output_filename))

    class AlignmentCheckerConfig:
        class NoGapsNoNAlignmentCheckerConfig:
            def __init__(self):
                self.gap_symbol = '-'

        def __init__(self):
            self.method = NoGapsNoNAlignmentChecker
            if self.method == NoGapsNoNAlignmentChecker:
                self.params = self.NoGapsNoNAlignmentCheckerConfig()

        def __str__(self):
            result = "AlignmentChecker: "
            if self.method == NoGapsNoNAlignmentChecker:
                result += "only sequences with no gaps are considered"
            result += "\n"
            return result

    class AlignmentCropperConfig:
        class UptoLastReliableKMerAlignmentCropperConfig:
            def __init__(self, command_line_args):
                self.k_mer_len = command_line_args.k_mer_length

        def __init__(self, command_line_args):
            self.method = UptoLastReliableKMerAlignmentCropper
            if self.method == UptoLastReliableKMerAlignmentCropper:
                self.params = self.UptoLastReliableKMerAlignmentCropperConfig(command_line_args)

        def __str__(self):
            result = "AlignmentCropper: "
            if self.method == UptoLastReliableKMerAlignmentCropper:
                result += "crop alignment up to last exact %d-mer match" % \
                          (self.params.k_mer_len)
            result += "\n"
            return result

    def __init__(self, command_line_args):
        self.IO = self.IOConfig(command_line_args)
        self.alignment_checker = self.AlignmentCheckerConfig()
        self.alignment_cropper = self.AlignmentCropperConfig(command_line_args)
        self.k_mer_len = command_line_args.k_mer_length

    def __str__(self):
        output_str = "========== CONFIG ==========\n"
        output_str += self.IO.__str__()
        output_str += self.alignment_checker.__str__()
        output_str += self.alignment_cropper.__str__()
        output_str += str(self.k_mer_len) + "-mers are used\n"
        output_str += "============================\n"
        return output_str
