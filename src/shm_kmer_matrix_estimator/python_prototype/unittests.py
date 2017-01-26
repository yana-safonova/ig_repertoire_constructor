#!/usr/bin/env python2

from __future__ import print_function

import unittest
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from alignment_checkers import NoGapsAlignmentChecker
from config import Config
from germline_alignment import GermlineAlignment
from mismatch_finder import NoKNeighboursMismatchFinder


class TestNoKNeighboursMismatchFinder(unittest.TestCase):
    config = Config.MismatchFinderConfig.NoKNeighboursMismatchFinderConfig(k_mer_len=5)
    finder = NoKNeighboursMismatchFinder(config)

    def test(self):
        self.assertEqual(self.finder.find_mismatch_positions('AACCGGTTAA',
                                                             'AATCGGAAAA'), [2])

        self.assertEqual(self.finder.find_mismatch_positions('AACTGGTTAA',
                                                             'AATCGGAAAA'), [])

        self.assertEqual(self.finder.find_mismatch_positions(
            'TCTCCAACGTTTTCTGTGCACGAGGGAGGT',
            'TCTCCAACGTTTTCTGTGCACGAGGGAGGT'), [])

        self.assertEqual(self.finder.find_mismatch_positions(
            'TCTCCAACGTTTTCTGTGCACGAGGGAGGT',
            'CCTCCATCGGTTTCTGTGCACGAGGGAGGT'), [6, 9])

        self.assertEqual(self.finder.find_mismatch_positions(
            'TTTCCAACGTTTTCTGTGCACGAGGGAGGT',
            'CCTCCATCGGTTTCTGTGCACGAGGGAGGT'), [6, 9])

        self.assertEqual(self.finder.find_mismatch_positions(
            'TTTCCAACGTTTTCTGTGCACGAGGGAGGT',
            'CCTCCATCGGTTTCTGTGCATGAGGGAGGT'), [6, 9, 20])

        self.assertEqual(self.finder.find_mismatch_positions(
            'TTTCCAACGTTTTCTGTGCACGAGGGAGGT',
            'CCTCCATCGGTTTCTGTGCATCAGGGAGGT'), [6, 9])

        self.assertEqual(self.finder.find_mismatch_positions(
            'TTTCCAACGTTTTCTGTGCACGAGGGAGGT',
            'CCTCCATCGGTTTCTGTGCATGACGGAGGT'), [6, 9, 20, 23])

        self.assertEqual(self.finder.find_mismatch_positions('CACCGGAAAA',
                                                             'AACTGGAAAA'), [3])

        self.assertEqual(self.finder.find_mismatch_positions('CACCGGAAAA',
                                                             'AAACGGAAAA'), [])


class TestNoGapsAlignmentChecker(unittest.TestCase):
    config = Config.AlignmentCheckerConfig.NoGapsAlignmentCheckerConfig()
    checker = NoGapsAlignmentChecker(config)

    def test(self):
        rec1 = SeqRecord(Seq('AGTACACTGGT', generic_dna))
        rec2 = SeqRecord(Seq('AGTAC-CTGGT', generic_dna))
        self.assertFalse(self.checker.check(GermlineAlignment(rec1, rec2)))

        rec1 = SeqRecord(Seq('AGTACACTGGT', generic_dna))
        rec2 = SeqRecord(Seq('AGTACCCTGGT', generic_dna))
        self.assertTrue(self.checker.check(GermlineAlignment(rec1, rec2)))

        rec1 = SeqRecord(Seq('AGTAC-CTGGT', generic_dna))
        rec2 = SeqRecord(Seq('AGTAC-CTGGT', generic_dna))
        self.assertFalse(self.checker.check(GermlineAlignment(rec1, rec2)))

        rec1 = SeqRecord(Seq('AGTAC-CTGGT', generic_dna))
        rec2 = SeqRecord(Seq('AGTACCCTGGT', generic_dna))
        self.assertFalse(self.checker.check(GermlineAlignment(rec1, rec2)))

        rec1 = SeqRecord(Seq('---', generic_dna))
        rec2 = SeqRecord(Seq('---', generic_dna))
        self.assertFalse(self.checker.check(GermlineAlignment(rec1, rec2)))


if __name__ == '__main__':
    unittest.main()
