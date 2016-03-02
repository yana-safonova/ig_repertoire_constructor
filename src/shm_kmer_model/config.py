#!/usr/bin/env python2

from __future__ import print_function

from germline_alignment import NoGapsAlignmentChecker

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

    def __init__(self, command_line_args):
        self.IO = self.IOConfig(command_line_args)
        self.alignment_checker = NoGapsAlignmentChecker #TODO from config.file

    def __str__(self):
        output_str = "========== CONFIG ==========\n"
        output_str += self.IO.__str__()
        output_str += "============================\n"
        return output_str
