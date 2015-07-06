#!/usr/bin/env python

import argparse
from collections import defaultdict
import mass_spectra_analysis.file_utils as file_utils
import os
import os.path
import sys
import itertools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

def DrawCoverageDistribution(data, filename, plotname):
    last_i = len(data)
    for i in range(len(data) - 1, 0, -1):
        if data[i] != 0:
            last_i = i + 1
            break
    data = data[:last_i]
    data.append(0)
    fig, ax = plt.subplots()
    # ind = np.arange(len(data))
    cdr_color = "gray"
    height = max(data) * 1.1
    plt.gca().add_patch(matplotlib.patches.Rectangle((27, 0), 12, height, facecolor=cdr_color, lw=0))
    plt.gca().add_patch(matplotlib.patches.Rectangle((56, 0), 10, height, facecolor=cdr_color, lw=0))
    plt.gca().add_patch(matplotlib.patches.Rectangle((95, 0), 20, height, facecolor=cdr_color, lw=0))
    plt.fill(range(len(data)), data, 'b')
    plt.xlabel('Position')
    plt.ylabel('Coverage')
    plt.title(plotname)
    plt.savefig(filename)
    plt.close()

def DrawHistogram(data, xlabel, ylabel, plotname, filename):
    if not data:
        return
    fig, ax = plt.subplots()
    ind = np.arange(len(data))
    n, bins, patches = plt.hist(data, 20, facecolor='b')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(plotname)
    plt.savefig(filename)
    plt.close()

class Metrics:
    metric_names = ['Aligned scans', 'PSMs', 'Peptides', 'Covered sequences', \
                    'Average PSM number per scan', 'Maximum PSM number per scan', \
                    'Average peptide length', 'Minimum peptide length', 'Maximum peptide length',
                    'Average PSM number per covered sequence', 'Maximum PSM number per covered sequence'] # Peptides - uniquely aligned scans

    def __init__(self, mass_spec_alns):
        self.mass_spec_alns = mass_spec_alns
        self.regions = None
        self.general_metrics = defaultdict(dict) # key - spectra name or 'Total'
        self.PSM_coverage_distr = {}
        self.peptide_coverage_distr = {}
        self.peptide_lengths = defaultdict(list)
        self.sequences_coverage = defaultdict(dict) # seq_id -> {spectra_name -> [identifications]}
        self.covered_cdrs = {} # spectra_name -> {reg_name -> sequences that have at least one PSM on this CDR}
        self.PSM_reg_coverage = {} # spectra_name -> {reg_name -> PSMs on this region}
        self.region_names = set()

    def InitMassSpectraStats(self, spectra_name):
        self.PSM_coverage_distr[spectra_name] = [0] * 300
        self.peptide_coverage_distr[spectra_name] = [0] * 300
        self.general_metrics[spectra_name]['Aligned scans'] = 0
        self.general_metrics[spectra_name]['PSMs'] = 0
        self.general_metrics[spectra_name]['Peptides'] = 0
        self.sequences_coverage[spectra_name] = {}
        if spectra_name != 'Total':
            self.PSM_reg_coverage[spectra_name] = {}
            self.covered_cdrs[spectra_name] = {}

    def GetAllSpectraNames(self):
        return (aln.spectra_name for aln in self.mass_spec_alns)

    def CalculateTotalMetrics(self):
        if len(self.mass_spec_alns) == 1:
            return
        self.general_metrics['Total']['Aligned scans'] = \
            sum(self.general_metrics[spectra_name]['Aligned scans']
                for spectra_name in self.GetAllSpectraNames())
        self.general_metrics['Total']['PSMs'] = \
            sum(self.general_metrics[spectra_name]['PSMs']
                for spectra_name in self.GetAllSpectraNames())
        self.general_metrics['Total']['Peptides'] = \
            sum(self.general_metrics[spectra_name]['Peptides']
                for spectra_name in self.GetAllSpectraNames())
        self.general_metrics['Total']['Covered sequences'] = \
            len(self.sequences_coverage)
        self.general_metrics['Total']['Average PSM number per scan'] = \
            sum(sum(len(v) for v in aln.spectrum_identifications.values())
                for aln in self.mass_spec_alns) / \
            float(sum(len(aln.spectrum_identifications)
                for aln in self.mass_spec_alns))
        self.general_metrics['Total']['Maximum PSM number per scan'] = \
            max(self.general_metrics[spectra_name]['Maximum PSM number per scan']
                for spectra_name in self.GetAllSpectraNames())
        self.general_metrics['Total']['Average peptide length'] = \
            sum(sum(v) for v in self.peptide_lengths.values()) / \
            float(sum(len(v) for v in self.peptide_lengths.values()))
        self.general_metrics['Total']['Minimum peptide length'] = \
            min(self.general_metrics[spectra_name]['Minimum peptide length'] \
                for spectra_name in self.GetAllSpectraNames())
        self.general_metrics['Total']['Maximum peptide length'] = \
            max(self.general_metrics[spectra_name]['Maximum peptide length'] \
                for spectra_name in self.GetAllSpectraNames())
        self.peptide_lengths['Total'] = list(itertools.chain(*self.peptide_lengths.values()))

    def CalculateSingleMetrics(self, mass_spec_aln):
        spectra_name = mass_spec_aln.spectra_name
        self.general_metrics[spectra_name]['Average PSM number per scan'] = \
            sum(len(v) for v in mass_spec_aln.spectrum_identifications.values()) /\
            float(len(mass_spec_aln.spectrum_identifications))
        self.general_metrics[spectra_name]['Maximum PSM number per scan'] = \
            max(len(v) for v in mass_spec_aln.spectrum_identifications.values())
        self.general_metrics[spectra_name]['Average peptide length'] = \
            sum(self.peptide_lengths[spectra_name]) / \
            float(len(self.peptide_lengths[spectra_name]))
        self.general_metrics[spectra_name]['Maximum peptide length'] = \
            max(self.peptide_lengths[spectra_name])
        self.general_metrics[spectra_name]['Minimum peptide length'] = \
            min(self.peptide_lengths[spectra_name])
        self.general_metrics[spectra_name]['Covered sequences'] = \
            sum((spectra_name in seq_identifications) for seq_identifications in \
                self.sequences_coverage.values())

    def ProcessSingleAlignment(self, mass_spec_aln):
        spectra_name = mass_spec_aln.spectra_name

        for spectrum_id, identifications in mass_spec_aln.spectrum_identifications.items():
            self.general_metrics[spectra_name]['Aligned scans'] += 1
            self.general_metrics[spectra_name]['PSMs'] += len(identifications)
            if len(identifications) == 1:
                self.general_metrics[spectra_name]['Peptides'] += 1
            for identification in identifications:
                peptide_length = identification.seq_to - identification.seq_from + 1
                self.peptide_lengths[spectra_name].append(peptide_length)
                if spectra_name not in self.sequences_coverage[identification.seq_id]:
                    self.sequences_coverage[identification.seq_id][spectra_name] = []
                self.sequences_coverage[identification.seq_id][spectra_name].append(identification)
                for pos in range(identification.seq_from, identification.seq_to + 1):
                    self.PSM_coverage_distr[spectra_name][pos] += 1
                    self.PSM_coverage_distr['Total'][pos] += 1
                    if len(identifications) == 1:
                        self.peptide_coverage_distr[spectra_name][pos] += 1
                        self.peptide_coverage_distr['Total'][pos] += 1

    def CalculatePSMPerSeqMetrics(self):
        psm_per_seq = defaultdict(list) # spectra_name -> [PSM counts]
        for seq_id, data in self.sequences_coverage.items():
            for spectra_name, identifications in data.items():
                psm_per_seq[spectra_name].append(len(identifications))

        for aln in self.mass_spec_alns:
            spectra_name = aln.spectra_name
            self.general_metrics[spectra_name]['Average PSM number per covered sequence'] = \
                float(sum(psm_per_seq[spectra_name])) / \
                len(psm_per_seq[spectra_name])
            self.general_metrics[spectra_name]['Maximum PSM number per covered sequence'] = max(psm_per_seq[spectra_name])
        self.general_metrics['Total']['Average PSM number per covered sequence'] =\
            float(sum([sum(v) for v in psm_per_seq.values()])) / \
            sum([len(v) for v in psm_per_seq.values()])
        self.general_metrics['Total']['Maximum PSM number per covered sequence'] =\
            max([max(v) for v in psm_per_seq.values()])

    def GetCDRForSpectraIdentification(self, curr_regions, identification):
        MATCH_THRESHOLD = 2
        for reg_name, (start, end) in curr_regions.items():
            if not reg_name.startswith('CDR'):
                continue
            peptide_length = identification.seq_to - identification.seq_from + 1
            match_length = min(identification.seq_to, end) - max(identification.seq_from, start) + 1
            if match_length >= peptide_length / 2 and \
                    match_length >= MATCH_THRESHOLD:
                return reg_name
        return 'FR'

    def CalculateRegionMetrics(self):
        if not self.regions:
            return
        self.PSM_reg_coverage['Any'] = {}
        self.PSM_reg_coverage['Any'] = {}
        self.covered_cdrs['Any'] = {}
        for seq_id, identifications_dict in self.sequences_coverage.items():
            if seq_id not in self.regions.regions_dict:
                continue
            curr_regions = self.regions.regions_dict[seq_id]
            for spectra_name, identifications in identifications_dict.items():
                seq_regions = set()
                for identification in identifications:
                    reg_name = self.GetCDRForSpectraIdentification(curr_regions, identification)
                    seq_regions.add(reg_name)
                    if reg_name not in self.PSM_reg_coverage[spectra_name]:
                        self.PSM_reg_coverage[spectra_name][reg_name] = 0
                    if reg_name not in self.PSM_reg_coverage['Any']:
                        self.PSM_reg_coverage['Any'][reg_name] = 0
                    self.PSM_reg_coverage[spectra_name][reg_name] += 1
                    self.PSM_reg_coverage['Any'][reg_name] += 1

                for reg_name in seq_regions:
                    self.region_names.add(reg_name)
                    if reg_name not in self.covered_cdrs[spectra_name]:
                        self.covered_cdrs[spectra_name][reg_name] = 0
                    if reg_name not in self.covered_cdrs['Any']:
                        self.covered_cdrs['Any'][reg_name] = 0
                    self.covered_cdrs[spectra_name][reg_name] += 1
                    self.covered_cdrs['Any'][reg_name] += 1
        self.region_names = sorted(self.region_names)

    def CalculateMetrics(self):
        self.InitMassSpectraStats('Total')
        for mass_spec_aln in self.mass_spec_alns:
            self.InitMassSpectraStats(mass_spec_aln.spectra_name)
            self.ProcessSingleAlignment(mass_spec_aln)
            self.CalculateSingleMetrics(mass_spec_aln)
        self.CalculatePSMPerSeqMetrics()
        self.CalculateTotalMetrics()

    def OutputMetrics(self, dirname):
        spectra_names = list(self.GetAllSpectraNames())
        print spectra_names

        handler = open(os.path.join(dirname, 'metrics.txt'), 'w')
        curr_names = [a for a in spectra_names]
        if len(spectra_names) != 1:
            curr_names.append('Total')
        handler.write('\t' + '\t'.join(curr_names) + '\n')
        for metric_name in self.metric_names:
            handler.write(metric_name + '\t' + \
                '\t'.join(str(self.general_metrics[spectra_name][metric_name])
                        if metric_name in self.general_metrics[spectra_name] else '*'
                    for spectra_name in curr_names) + '\n')
        handler.close()

        handler = open(os.path.join(dirname, 'covered_cdrs.txt'), 'w')
        handler.write('\t' + '\t'.join(self.region_names) + '\n')
        for spectra_name in spectra_names + ['Any']:
            if spectra_name == 'Any' and len(self.mass_spec_alns) == 1:
                continue
            data = self.covered_cdrs[spectra_name]
            handler.write(spectra_name + '\t' + \
                '\t'.join(str(data[reg_name]) if reg_name in data else '0'
                    for reg_name in self.region_names) + '\n')
        handler.close()

        handler = open(os.path.join(dirname, 'psm_on_ig_regions.txt'), 'w')
        handler.write('\t' + '\t'.join(self.region_names) + '\n')
        for spectra_name in spectra_names + ['Any']:
            data = self.PSM_reg_coverage[spectra_name]
            if spectra_name == 'Any' and len(self.mass_spec_alns) == 1:
                continue
            handler.write(spectra_name + '\t' + \
                '\t'.join(str(data[reg_name]) if reg_name in data else '0'
                    for reg_name in self.region_names) + '\n')
        handler.close()

        for spectra_name in spectra_names + ['Total']:
            distr = self.PSM_coverage_distr[spectra_name]
            if spectra_name == 'Total' and len(self.mass_spec_alns) == 1:
                continue
            filename = os.path.join(dirname, 'PSM_cov_' + os.path.basename(spectra_name) + '.png')
            DrawCoverageDistribution(distr, filename, 'PSM coverage')

        for spectra_name in spectra_names + ['Total']:
            distr = self.peptide_coverage_distr[spectra_name]
            if spectra_name == 'Total' and len(self.mass_spec_alns) == 1:
                continue
            filename = os.path.join(dirname, 'peptide_cov_' + os.path.basename(spectra_name) + '.png')
            DrawCoverageDistribution(distr, filename, 'Peptide coverage')

        for mass_spec_aln in self.mass_spec_alns:
            filename = os.path.join(dirname, 'PSM_per_scan_' + os.path.basename(mass_spec_aln.filename) + '.png')
            data = [len(v) for v in mass_spec_aln.spectrum_identifications.values()]
            threshold = np.percentile(data, [98])[0]
            DrawHistogram([v for v in data if v < threshold], \
                'PSMs per scan', 'Count', 'PSM number per scan', filename)

        for spectra_name in spectra_names + ['Total']:
            if spectra_name == 'Total' and len(self.mass_spec_alns) == 1:
                continue
            length_list = self.peptide_lengths[spectra_name]
            if spectra_name == 'Total' and len(self.mass_spec_alns) == 1:
                continue
            filename = os.path.join(dirname, 'peptide_length_' + os.path.basename(spectra_name) + '.png')
            DrawHistogram(length_list, 'Peptide length', 'Count', 'Peptide length distribution', filename)

def PrepareArguments():
    home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    parser = argparse.ArgumentParser(description = 'mass spectra analysis')
    parser.add_argument('--regions',
                        help='regions file')
    parser.add_argument('--outdir',
                        dest='output_dir',
                        help='output directory')
    parser.add_argument('inputs',
                        help='mass spectrum alignment files in mzIdentML 1.1 format',
                        nargs='*')
    parser.add_argument('--test',
                        help='run on the test dataset',
                        action='store_true')
    args = parser.parse_args()
    if args.test:
        args.output_dir = os.path.join(home_directory, 'test_output')
        args.regions = os.path.join(home_directory,
                                    'mass_spectra_analysis',
                                    'test_dataset',
                                    'example_regions.txt')
        args.inputs = [os.path.join(home_directory,
                                    'mass_spectra_analysis',
                                    'test_dataset',
                                    'example_HC_chymo_CID.mzid.spectra'),
                       os.path.join(home_directory,
                                    'mass_spectra_analysis',
                                    'test_dataset',
                                    'example_HC_trypsin_CID.mzid.spectra')]
    if not args.inputs:
        print 'You must specify at least one input file'
        parser.print_help()
        sys.exit(-1)
    if not args.output_dir:
        print 'You must specify output directory'
        parser.print_help()
        sys.exit(-1)
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    return args

if __name__ == '__main__':
    args = PrepareArguments()
    mass_spec_alns = []

    for filename in args.inputs:
        print 'Parsing ' + filename
        aln = file_utils.MassSpectraAlignment()
        if filename.endswith('.spectra'):
            aln.ParseSpectra(filename)
        else:
            aln.ParseMzindentFile(filename)
        mass_spec_alns.append(aln)
        print 'Parsing ' + filename + ' done'

    metrics = Metrics(mass_spec_alns)
    if args.regions:
        print 'Parsing regions file'
        regions = file_utils.Regions()
        regions.ParseIgBlastAlignment(args.regions)
        metrics.regions = regions
        print 'Parsing regions file done'

    metrics.CalculateMetrics()
    metrics.CalculateRegionMetrics()
    metrics.OutputMetrics(args.output_dir)
