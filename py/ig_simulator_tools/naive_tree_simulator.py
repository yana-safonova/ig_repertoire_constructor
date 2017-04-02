import os
from Bio import SeqIO
import numpy as np
import random

current_dir = os.path.dirname(__file__)
root_dir = os.path.realpath(os.path.join(current_dir, "../../"))
input_file = os.path.join(root_dir, "ig_simulator_test/test.fa")

def random_mutations(antibody, pois_lamb=3):
    antibody = list(antibody)
    n_mut = np.random.poisson(pois_lamb, 1)[0]
    mut_ind = np.random.choice(range(len(antibody)), size=n_mut, replace=False)
    for i in mut_ind:
        bases = list("ACGT")
        bases.remove(antibody[i])
        antibody[i] = bases[np.random.randint(3)]
    return "".join(antibody)


class Forest(object):
    def __init__(self, record):
        record.name = str(record.name)
        self.name, self.multiplicity = [x for x in record.name.split('_') if x.isdigit()]
        self.metaroot_seq = str(record.seq)


class Node(object):
    def __init__(self, seq, numb):
        self.seq = seq
        self.children = []
        self.numb = numb


def generate_tree(seq, n, ret_prob=0.5):
    pool = []
    root = Node(seq, 0)
    pool.append(root)
    for i in xrange(n):
        e = random.choice(pool)
        mut_e = Node(random_mutations(e.seq), i + 1)
        e.children.append(mut_e)
        bern = np.random.binomial(1, ret_prob, 1)[0]
        if not bern:
            pool.remove(e)
        pool.append(mut_e)
    return {'root': root, 'pool': pool}


def print_tree(root):
    if not root.children:
        return
    print(root.numb)
    for x in root.children:
        print(root.numb, x.numb)
        print_tree(x)


def go(forests, tree_size=100):
    results = [0] * len(forests)
    for i, forest in enumerate(forests):
        results[i] = []
        for _ in xrange(len(forest.multiplicity)):
            results[i].append(generate_tree(forest.metaroot_seq, tree_size))
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print_tree(results[i][-1]["root"])
    return results


def main():
    # np.random.seed(1)
    with open(input_file, "r") as f:
        records = list(SeqIO.parse(f, "fasta"))
    forests = [Forest(record) for record in records]
    results = go(forests)
    # root, pool = generate_tree(forests[0].metaroot_seq, 10000)
    # # print_tree(root)
    # print(len(pool))
    # seq_pool = [root.seq for root in pool]
    # print(seq_pool[-2])

if __name__ == "__main__":
    main()
