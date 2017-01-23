#!/usr/bin/env python2

from __future__ import print_function

from abc import ABCMeta, abstractmethod


class AbstractAlignmentChecker:
    __metaclass__ = ABCMeta

    @abstractmethod
    def check(self, alignment):
        pass


class NoGapsAlignmentChecker(AbstractAlignmentChecker):
    def check(self, alignment):
        return alignment.read.count(self.gap_symbol) == 0 and \
               alignment.germline_seq.count(self.gap_symbol) == 0

    def __str__(self):
        return "Alignment Checker: only sequences with no gaps are considered\n"

    def __init__(self, config):
        self.gap_symbol = config.gap_symbol
