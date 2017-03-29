import os.path

class Regions:
    def __init__(self):
        self.regions_dict = {} # seq_name -> {region_name -> (from, to)}

    def ParseIgBlastAlignment(self, filename):
        query = None
        regions = {}
        for l in open(filename):
            if l.startswith('Query='):
                if query:
                    self.regions_dict[query] = regions
                    regions = {}
                query = l.strip().split()[-1]
            if l.startswith('FR') or l.startswith('CDR'):
                fields = l.strip().split()
                regions[fields[0]] = (int(fields[1]), int(fields[2]))

    '''
    def ParseRegionsFile(self, filename):
        handler = open(filename)
        header_fields = handler.readline().strip().split()
        for l in handler:
            fields = l.strip().split()
            query = fields[0]
            from_, to = 0, 0
            self.regions_dict[query] = {}
            for i, value in enumerate(fields):
                if value == '*':
                    value = -1
                if header_fields[i].endswith('_from'):
                    from_ = int(value)
                elif header_fields[i].endswith('_to'):
                    to = int(value)
                    reg_name = '_'.join(header_fields[i].split('_')[:-1])
                    self.regions_dict[query][reg_name] = (from_, to)

    def SaveRegions(self, filename):
        handler = open(filename, 'w')
        for query, regions in self.regions_dict.items():
            handler.write('Query= ' + query + '\n')
            for region_name, (reg_from, reg_to) in regions.items():
                handler.write(region_name + '\t' + str(reg_from) + '\t' + str(reg_to) + '\n')
        handler.close()
    '''


class PeptideEvidence:
    def __init__(self, seq_id, start, end):
        self.seq_id = seq_id
        self.start = start
        self.end = end

class SeqSpectrumIdentification:
    def __init__(self, spectrum_id, start, end):
        self.spectrum_id = spectrum_id
        self.start = start
        self.end = end

class SpectrumIdentification:
    def __init__(self, seq_id, seq_from, seq_to):
        self.seq_id = seq_id
        self.seq_from = seq_from
        self.seq_to = seq_to

class MassSpectraAlignment:
    def __init__(self):
        self.filename = None
        self.spectra_name = None
        self.spectrum_identifications = {}
        self.regions = Regions()

    def ParseMzindentFile(self, filename): # by hand, may be should use pyteomics
        self.filename = filename
        self.spectra_name = os.path.basename(filename)
        dbseq_ids = {} # map from database IDs to original sequence IDs
        peptide_evid_ids = {}
        tree = ET.parse(filename)
        root = tree.getroot()
        xmlns = '{http://psidev.info/psi/pi/mzIdentML/1.1}'
                
        for dbseq in root.iter(xmlns + 'DBSequence'):
            dbseq_ids[dbseq.get('id')] = dbseq.get('accession')

        for peptide_evid in root.iter(xmlns + 'PeptideEvidence'):
            peptide_evid_ids[peptide_evid.get('id')] = PeptideEvidence(
                dbseq_ids[peptide_evid.get('dBSequence_ref')],
                peptide_evid.get('start'),
                peptide_evid.get('end'))

        seq_peptides = defaultdict(list)

        for result in root.iter(xmlns + 'SpectrumIdentificationResult'):
            spectrum_id = int(result.get('spectrumID').split('=')[1])
            identifications = []
            for item in result.findall(xmlns + 'SpectrumIdentificationItem'):
                for ref in item.findall(xmlns + 'PeptideEvidenceRef'):
                    pep_evid = peptide_evid_ids[ref.get('peptideEvidence_ref')]
                    identifications.append(SpectrumIdentification(pep_evid.seq_id, int(pep_evid.start), int(pep_evid.end)))
            self.spectrum_identifications[spectrum_id] = identifications

    def ParseSpectra(self, filename):
        self.filename = filename
        self.spectra_name = os.path.basename(filename)
        spectrum_id = None
        identifications = []
        for l in open(filename):
            fields = l.strip().split()
            if len(fields) == 2 and l.startswith('ID:'):
                if spectrum_id and identifications:
                    self.spectrum_identifications[spectrum_id] = identifications
                    identifications = []
                spectrum_id = int(l.strip().split()[1])
            else:
                identifications.append(SpectrumIdentification(fields[0], int(fields[1]), int(fields[2])))
        self.spectrum_identifications[spectrum_id] = identifications

    def SaveSpectra(self, filename):
        handler = open(filename, 'w')
        for spectrum_id, identifications in self.spectrum_identifications.items():
            handler.write('ID:\t' + str(spectrum_id) + '\n')
            for identification in identifications:
                handler.write('\t'.join((identification.seq_id, identification.seq_from, identification.seq_to)) + '\n')
        handler.close()


