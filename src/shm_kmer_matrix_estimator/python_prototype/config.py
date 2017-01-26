#!/usr/bin/env python2

from __future__ import print_function

from alignment_checkers import NoGapsAlignmentChecker
from alignment_cropper import UptoLastReliableKMerAlignmentCropper
from mismatch_finder import NoKNeighboursMismatchFinder
from mismatch_finder import TrivialMismatchFinder


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
        class NoGapsAlignmentCheckerConfig:
            def __init__(self):
                self.gap_symbol = '-'

        def __init__(self):
            self.method = NoGapsAlignmentChecker
            if self.method == NoGapsAlignmentChecker:
                self.params = self.NoGapsAlignmentCheckerConfig()

        def __str__(self):
            result = 'AlignmentChecker: '
            if self.method == NoGapsAlignmentChecker:
                result += 'only sequences with no gaps are considered'
            result += '\n'
            return result

    class AlignmentCropperConfig:
        class UptoLastReliableKMerAlignmentCropperConfig:
            def __init__(self, command_line_args):
                self.k_mer_len = command_line_args.k_mer_len

        def __init__(self, command_line_args):
            self.method = UptoLastReliableKMerAlignmentCropper
            if self.method == UptoLastReliableKMerAlignmentCropper:
                self.params = self.UptoLastReliableKMerAlignmentCropperConfig(command_line_args)

        def __str__(self):
            result = 'AlignmentCropper: '
            if self.method == UptoLastReliableKMerAlignmentCropper:
                result += 'crop alignment up to last exact %d-mer match' % \
                          (self.params.k_mer_len)
            result += '\n'
            return result

    class MismatchFinderConfig:
        class NoKNeighboursMismatchFinderConfig:
            def __init__(self, params=None, k_mer_len=None):
                if params is not None:
                    self.k_mer_len = params.k_mer_len
                if k_mer_len is not None:
                    self.k_mer_len = k_mer_len

        def __init__(self, command_line_args):
            if command_line_args.mismatch_finder == 'trivial':
                self.method = TrivialMismatchFinder
                self.params = None
            elif command_line_args.mismatch_finder == 'k_neighbour':
                self.method = NoKNeighboursMismatchFinder
                self.params = self.NoKNeighboursMismatchFinderConfig(command_line_args)

    def __init__(self, command_line_args):
        self.IO = self.IOConfig(command_line_args)
        self.alignment_checker = self.AlignmentCheckerConfig()
        self.alignment_cropper = self.AlignmentCropperConfig(command_line_args)
        self.mismatch_finder = self.MismatchFinderConfig(command_line_args)
        self.k_mer_len = command_line_args.k_mer_len

    def __str__(self):
        output_str = '========== CONFIG ==========\n'
        output_str += self.IO.__str__()
        output_str += self.alignment_checker.__str__()
        output_str += self.alignment_cropper.__str__()
        output_str += str(self.k_mer_len) + '-mers are used\n'
        output_str += '============================\n'
        return output_str
