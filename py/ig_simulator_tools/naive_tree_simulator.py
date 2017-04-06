from __future__ import division
import os
import errno

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import random


def smart_makedirs(dirname):
    try:
        os.makedirs(dirname)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise exc


def random_mutations(antibody, pois_lamb=3):
    antibody = list(antibody)
    n_mut = np.random.poisson(pois_lamb, 1)[0]
    mut_ind = np.random.choice(range(len(antibody)), size=n_mut, replace=False)
    for i in mut_ind:
        if antibody[i] == "N":
            continue
        bases = list("ACGT")
        bases.remove(antibody[i])
        antibody[i] = bases[np.random.randint(3)]
    return "".join(antibody)


class ParsedRecord(object):
    def __init__(self, record):
        record.name = str(record.name)
        self.name, self.multiplicity = [int(x) for x in record.name.split('_') if x.isdigit()]
        self.metaroot_seq = str(record.seq)


class Node(object):
    def __init__(self, seq, numb):
        self.seq = seq
        self.children = []
        self.numb = numb
        self.included = True


def generate_tree(seq, n, ret_prob):
    pool = []
    root = Node(seq, 0)
    pool.append(root)
    indeces = []
    weights = []
    for i in xrange(n):
        indeces.append(i)
        weights.append(1)
        index = np.random.choice(indeces, size=1,
                                 p=np.array(weights) / sum(weights))[0]
        weights[index] += 1
        e = pool[index]
        mut_e = Node(random_mutations(e.seq), i + 1)
        e.children.append(mut_e)
        bern = np.random.binomial(1, ret_prob, 1)[0]
        if not bern:
            e.included = False
        pool.append(mut_e)

    results = {'root': root}
    results['all_seqs'] = [(id, x.seq) for id, x in enumerate(pool)]
    results['filtered_seqs'] = [(id, x.seq) for id, x in enumerate(pool) if x.included]
    results['edge_list'] = edge_list(root)
    return results


def edge_list(root):
    elist = [(root.numb, x.numb) for x in root.children]
    for x in root.children:
        elist += edge_list(x)
    return elist


def go(records, lamb, ret_prob):
    records = [ParsedRecord(record) for record in records]
    results = []
    for i, record in enumerate(records):
        print(i + 1, len(records))
        results.append([generate_tree(record.metaroot_seq,
                                      np.random.geometric(lamb, size=1)[0],
                                      ret_prob=ret_prob)
                        for _ in xrange(record.multiplicity)])
    return results


def output_forests(forests, output_folder = ""):
    all_records, filtered_records = [], []

    def create_records(x):
        return [SeqRecord(Seq(seq),
                          id="metaroot_%d_tree_%d_id_%d" % (i, m, id),
                          description="")
                for id, seq in x]

    for i, forest in enumerate(forests):
        for m, tree in enumerate(forest):
            all_records += create_records(tree['all_seqs'])
            filtered_records += create_records(tree['filtered_seqs'])

    smart_makedirs(output_folder)
    SeqIO.write(filtered_records, os.path.join(output_folder, "filtered_records.fasta"), "fasta")
    SeqIO.write(all_records, os.path.join(output_folder, "all_records.fasta"), "fasta")

    edge_lists_dir = os.path.join(output_folder, "trees_edge_lists")
    smart_makedirs(edge_lists_dir)
    for i, forest in enumerate(forests):
        for m, tree in enumerate(forest):
            with open(os.path.join(edge_lists_dir, "antibody_%d_tree_%d.txt" %(i, m)), "w") as f:
                for tup in tree['edge_list']:
                    f.write("%d %d\n" % tup)


def ParseCommandLineParams():
    import argparse
    current_dir = os.path.dirname(os.path.realpath(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outdir",
                        type=str,
                        default=current_dir,
                        help="Output directory")
    root_dir = os.path.realpath(os.path.join(current_dir, "../../"))
    input_file = os.path.join(root_dir, "ig_simulator_test/test.fa")
    parser.add_argument("-i", "--input",
                        type=str,
                        default=input_file)
    parser.add_argument("--seed",
                        type=int,
                        default=int(np.random.randint(low=0, high=100000, size=1)[0]))
    parser.add_argument("-s", "--exp_tree_size",
                        type=float,
                        default=0.01,
                        help="Mean tree size")
    parser.add_argument("-p", "--ret_prob",
                        type=float,
                        default=0.5,
                        help="Probability to return chosen antibody to the pool")
    return parser.parse_args()


def main():
    params = ParseCommandLineParams()
    np.random.seed(params.seed)
    with open(params.input, "r") as f:
        records = list(SeqIO.parse(f, "fasta"))

    results = go(records, lamb=params.exp_tree_size, ret_prob=params.ret_prob)
    output_forests(results, params.outdir)


if __name__ == "__main__":
    main()
