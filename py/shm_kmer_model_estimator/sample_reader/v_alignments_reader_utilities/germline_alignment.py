#!/usr/bin/env python2

from __future__ import print_function


class GermlineAlignment(object):
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
