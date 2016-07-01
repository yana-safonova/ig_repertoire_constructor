import os
from Bio import SeqIO
import random
import numpy as np
from Bio import Seq
from rand import rcm2Rand


def mult2mult(clustering_filename, reference_filename, match_filename):
    from rep_sn import parse_multiplicity, parse_size

    with open(clustering_filename) as f:
        clustering_mults = [parse_size(rec.id) for rec in SeqIO.parse(f, "fasta")]

    with open(reference_filename) as f:
       reference_mults = [parse_multiplicity(rec.id) for rec in SeqIO.parse(f, "fasta")]

    unmatched_references = set(range(len(reference_mults)))
    result = []

    with open(match_filename) as f:
        for i, line in enumerate(f):
            line = line.strip()
            cl_mult = clustering_mults[i]
            ref_mult = 0
            if line:
                dist = -int(line.split()[0])
                if dist == 0:
                    neibs = map(int, line.split()[1:])
                    if len(neibs):
                        ref_mult = sum(reference_mults[j] for j in neibs)
                        unmatched_references.difference_update(neibs)

            result.append((cl_mult, ref_mult))

    for i in unmatched_references:
        result.append((0, reference_mults[i]))

    return result


def match_stat(fname, tau=2):
    table = [0] * (tau+1)
    N = 0
    with open(fname) as fh:
        for _i, line in enumerate(fh):
            if line.strip():
                dist = -int(line.strip().split()[0])
                if dist <= tau:
                    table[dist] += 1
            ###### else:
            ######     print _i

            N += 1
    return table, N


def pretty_stat(fname1, fname2, f, tau=4):
    tbl1, N1 = match_stat(fname1, tau)
    tbl2, N2 = match_stat(fname2, tau)

    res = []

    res.append(N1)
    print >>f,"Identified Igs: (from %d)" % N1
    s = 0
    for i in range(tau+1):
        s += tbl1[i]
        print >>f,"With distance %d: %d" % (i, s)
        res.append(s)

    res.append(N2)
    print >>f, "Extraneous sequences: (from %d)" % N2
    s = 0
    for i in range(tau+1):
        s += tbl2[i]
        print >>f, "With distance %d: %d" % (i, N2 - s)
        res.append(N2-s)

    res = ",".join(map(str, res))

    return res



def convert_mixcr_igrc(input_file, output_file):
    with open(input_file) as fh, open(output_file, "w") as fout:
        # Skip header
        fh.next()

        for i, line in enumerate(fh):
            seq, size = line.strip().split()
            size = int(size)
            fout.writelines(">cluster___%d___size___%d\n" % (i + 1, size))
            fout.writelines(seq)
            fout.writelines("\n")

def fa2fq(input_file, output_file):
    with open(input_file) as fh, open(output_file, "w") as fout:
        for record in SeqIO.parse(fh, "fasta"):
            # record.letter_annotations["phred_quality"] = [40] * len(record)
            record.letter_annotations["phred_quality"] = [random.randint(30, 50) for _ in xrange(len(record))]
            SeqIO.write(record, fout, "fastq")


def RC(l):
    S = set(list("ACTG"))
    s = S.difference([l])
    return random.choice(list(s))

def simulate(input_file, output_file, coverage=25, errors=1, maxsn=25, random_errors=True):
    with open(input_file) as fh, open(output_file, "w") as fout:
        for record in SeqIO.parse(fh, "fasta"):
            SS = str(record.seq)[:]
            desc = str(record.description).replace("abundance:", "")
            _maxsn = maxsn
            for i in xrange(coverage):
                rec = record
                n_errors = np.random.poisson(errors, 1)[0] if random_errors else errors

                if n_errors == 0:
                    if _maxsn:
                        _maxsn -= 1
                    else:
                        n_errors += 1

                positions = random.sample(range(min(len(rec.seq), 300)), n_errors)

                s = SS[:]
                s = list(s)
                for pos in positions:
                    s[pos] = RC(s[pos])
                rec.letter_annotations = {}
                rec.seq = Seq.Seq("".join(s))

                rec.description = desc + "__%d" % i
                rec.id = desc + "__%d" % i
                rec.letter_annotations["phred_quality"] = [random.randint(30, 50) for _ in xrange(len(rec))]

                SeqIO.write(rec, fout, "fastq")

# input_file = "/ssd/ig_repertoire_constructor/test_dataset/test7.fastq"
loci = "IGH"

mixcr_command = "java -jar mixcr.jar "

igrc_path = "/home/ashlemov/Git/ig_repertoire_constructor/"

igsim_path = "/home/ashlemov/Git/ig_simulator/"

trie_comp = "%s/build/release/bin/ig_trie_compressor " % igrc_path
matcher = "%s/build/release/bin/ig_matcher " % igrc_path



def F(coverage=25, errors=1, maxsn=25, out_csv="out.csv", random_errors=True):
    # fa2fq("final_repertoire.fa", "final_repertoire.fq")
    input_file = "input.fq"

    os.system("%s -i final_repertoire.fa -o reference.fa" % trie_comp)
    # simulate("reference.fa", input_file, coverage, errors, maxsn)
    simulate("final_repertoire.fa", input_file, coverage, errors, maxsn)

    # os.system("%s/igrec.py --debug --tau=2 -t 36 --loci IGH -s %s -o igrc_out" % (igrc_path, input_file))
    os.system("%s/igrec.py --debug --tau=4 -t 36 --loci IGH -s %s -o igrc_out" % (igrc_path, input_file))

    # fa2fq("igrc_out/vj_finder/cleaned_reads.fa", "igrc_out/vj_finder/cleaned_reads.fq")

    os.system("%s align -f -g -r align_report.txt --loci %s --noMerge -OvParameters.geneFeatureToAlign=VTranscript %s mixcr.vdjca" % (mixcr_command, loci, input_file))
    os.system("%s assemble -f -r assemble_report.txt -OassemblingFeatures=\"{FR1Begin:CDR3End}\" mixcr.vdjca mixcr.clns" % mixcr_command)
    # os.system("%s/mixcr assemble -f -r assemble_report.txt -OassemblingFeatures=\"{FR1Begin:FR4End}\" mixcr.vdjca mixcr.clns" % mixcr_path)
    # os.system("%s/mixcr assemble mixcr.vdjca mixcr.clns -f" % mixcr_path)
    os.system("%s exportClones -pf preset.pf -f --no-spaces mixcr.clns mixcr.txt" % mixcr_command)
    convert_mixcr_igrc("mixcr.txt", "mixcr_final.fa")
    os.system("./rep_sn.py %s %s --limit=%d" % ("mixcr_final.fa", "mixcr_large.fa", 5))
    os.system("./rep_sn.py %s %s --limit=%d" % ("mixcr_final.fa", "mixcr_large3.fa", 3))
    os.system("./rep_sn.py %s %s --limit=%d" % ("mixcr_final.fa", "mixcr_large2.fa", 2))



    # compress outputs
    os.system("%s -i mixcr_final.fa -o mixcr_res.fa" % trie_comp)
    os.system("%s -i igrc_out/final_repertoire.fa -o igrc_res.fa" % trie_comp)
    os.system("%s -i igrc_out/final_repertoire_large.fa -o igrc_res_large.fa" % trie_comp)
    os.system("%s -i igrc_out/super_reads.fa -o igrc_res_super_reads.fa" % trie_comp)
    os.system("%s -i mixcr_large.fa -o mixcr_res_large.fa" % trie_comp)
    os.system("%s -i mixcr_large3.fa -o mixcr_res_large3.fa" % trie_comp)
    os.system("%s -i mixcr_large2.fa -o mixcr_res_large2.fa" % trie_comp)


    os.system("%s -i reference.fa -I mixcr_res.fa -o mixcr1.match -O mixcr2.match -k 10 --tau 4" % matcher)
    os.system("%s -i reference.fa -I mixcr_res_large.fa -o mixcr1large.match -O mixcr2large.match -k 10 --tau 4" % matcher)
    os.system("%s -i reference.fa -I mixcr_res_large3.fa -o mixcr1large3.match -O mixcr2large3.match -k 10 --tau 4" % matcher)
    os.system("%s -i reference.fa -I mixcr_res_large2.fa -o mixcr1large2.match -O mixcr2large2.match -k 10 --tau 4" % matcher)
    os.system("%s -i reference.fa -I igrc_res.fa -o igrc1.match -O igrc2.match -k 10 --tau 4" % matcher)
    os.system("%s -i reference.fa -I igrc_res_large.fa -o igrc1large.match -O igrc2large.match -k 10 --tau 4" % matcher)
    os.system("%s -i reference.fa -I igrc_res_super_reads.fa -o igrc1super.match -O igrc2super.match -k 10 --tau 4" % matcher)


    m2m_igrec = mult2mult("igrc_res.fa", "reference.fa", "igrc2.match")

    with open("m2m_igrec_snodes_%d_lam_%f.txt" % (supernodes, error_rate), "w") as f:
        for mc, mr in m2m_igrec:
            print >>f, "%d,%d" % (mc, mr)

    m2m_mixcr = mult2mult("mixcr_res.fa", "reference.fa", "mixcr2.match")
    with open("m2m_mixcr_snodes_%d_lam_%f.txt" % (supernodes, error_rate), "w") as f:
        for mc, mr in m2m_mixcr:
            print >>f, "%d,%d" % (mc, mr)

    with open("results.log", "a") as f, open(out_csv, "w") as csv:
        print >>csv, "Tool,Nref,i0,i1,i2,i3,i4,Nout,e0,e1,e2,e3,e4"

        print >>f, "Coverage: %d, supernodes: %d, error_rate: %f " % (coverage, maxsn, errors)
        print >>f, "MiXCR:"
        s = pretty_stat("mixcr1.match", "mixcr2.match", f)
        print >>csv, "%s,%s" % ("MiXCR", s)
        print >>f, "\n"

        print >>f, "MiXCR large (>=5):"
        s = pretty_stat("mixcr1large.match", "mixcr2large.match", f)
        print >>csv, "%s,%s" % ("MiXCR5", s)
        print >>f, "\n"

        print >>f, "MiXCR large (>=3):"
        s = pretty_stat("mixcr1large3.match", "mixcr2large3.match", f)
        print >>csv, "%s,%s" % ("MiXCR3", s)
        print >>f, "\n"

        print >>f, "MiXCR large (>=2):"
        s = pretty_stat("mixcr1large2.match", "mixcr2large2.match", f)
        print >>csv, "%s,%s" % ("MiXCR2", s)
        print >>f, "\n"

        print >>f, "IGRC:"
        s = pretty_stat("igrc1.match", "igrc2.match", f)
        print >>csv, "%s,%s" % ("IGRC", s)
        print >>f, "\n"

        print >>f, "IGRC large (>=5):"
        s = pretty_stat("igrc1large.match", "igrc2large.match", f)
        print >>csv, "%s,%s" % ("IGRC5", s)
        print >>f, "\n"

        # print >>f, "IGRC super reads (>=5):"
        # pretty_stat("igrc1super.match", "igrc2super.match", f)
        # print >>f, "\n============================================\n"

if __name__ == "__main__":
    random.seed(42)


    with open("results.log", "w"):
        pass
    with open("clustering_measures.txt", "w") as f:
        print >>f,"tool,supernodes,lambda,FM,rand,rand_adj,Jaccard,FM_large,rand_large,rand_adj_large,jaccard_large"

    # os.system("%s/ig_simulator.py --chain-type HC --num-bases 100 --num-mutated 1000 --repertoire-size 5000 -o sim --skip-drawing" % igsim_path)
    os.system("%s/ig_simulator.py --chain-type HC --num-bases 1000 --num-mutated 10000 --repertoire-size 50000 -o sim --skip-drawing" % igsim_path)

    os.system("cp %s/sim/merged_reads.fastq ./input_merged_reads.fastq" % igsim_path)
    os.system("cp %s/sim/final_repertoire.fasta ./final_repertoire.fa" % igsim_path)

    print ">>> Reads simulated!"



    for coverage in [1]:
        # for supernodes in [1000000]:
        #     for error_rate in [0.1]:
        #         csv_file = "results/results_cov_%d_super_%d_lam_%f.csv" % (coverage,
        #                                                                    supernodes,
        #                                                                    error_rate)
        #         F(coverage, error_rate, supernodes, csv_file)
        # for supernodes in [0, 1000000]:
        #     for error_rate in [0.1, 0.5, 1, 1.1, 1.25, 1.4, 1.5, 2, 2.5, 3, 3.5]:
        for supernodes in [1000000, 0]:
            for error_rate in [0.1, 0.5, 1, 1.1, 1.25, 1.4, 1.5, 2, 2.5, 3, 3.5]:
                csv_file = "results/results_cov_%d_super_%d_lam_%f.csv" % (coverage,
                                                                           supernodes,
                                                                           error_rate)
                F(coverage, error_rate, supernodes, csv_file)
                with open("clustering_measures.txt", "a") as cm:
                    print >>cm, "%s,%d,%f,%s" % ("IGREC", supernodes, error_rate, ",".join(map(str, rcm2Rand())))
        # F(10, 2, 10001000000)

