import argparse
from collections import defaultdict
from copy import copy
import numpy as np
from Bio import SeqIO

import sys
sys.path.append('..')
from kmer_utilities import kmer_utilities as kmer_utils


class EM:
    def crop(self, read, gene):
        assert len(read) == len(gene)
        i = 0
        j = len(read)
        while read[i:(i+self.k)] != gene[i:(i+self.k)]:
            i += 1
        while read[(j-self.k):j] != gene[(j-self.k):j]:
            j -= 1
        return read[i:j], gene[i:j]

    def read_pinp(self, input_file):
        pinp = []
        it = SeqIO.parse(input_file, "fasta")
        try:
            while True:
                read, gene = str(it.next().seq), str(it.next().seq)
                if '-' in read + gene:
                    continue
                read, gene = self.crop(read, gene)
                pinp.append((read, gene))
        except StopIteration:
            pass
        return pinp

    def read_input(self, paired_input):
        self.nmut, self.mut = [0] * self.all_k, defaultdict(int)
        self.init_matr = np.zeros((self.all_k, 4))
        for read, gene in paired_input:
            for i in xrange(len(read) - self.k):
                rkmer, gkmer = read[i:(i+5)], gene[i:(i+5)]
                if 'N' in gkmer:
                    continue
                ind = kmer_utils.kmer_index(gkmer)
                if gkmer[2] == rkmer[2]:
                    self.nmut[ind] += 1
                    self.init_matr[ind, 0] += 1
                elif 'N' not in rkmer:
                    self.mut[(rkmer, gkmer)] += 1
                    read_ind = list("ACGT").index(rkmer[2])
                    gene_ind = list("ACGT").index(gkmer[2])
                    matr_ind = read_ind if read_ind > gene_ind else read_ind + 1
                    self.init_matr[ind, matr_ind] += 1
        self.nmut = np.array(self.nmut)
    
    def e_step(self):
        self.q = defaultdict(int)
        for pair in self.mut:
            read, gene = pair
            kmers, inds = self.get_valid_kmers(pair)
            log_q = np.log(self.nu[inds]) + np.log(self.theta[inds]) + \
                    np.log(self.pkb[inds, self.get_pkb_ind(inds, read[2])])
            self.q[pair] = np.exp(log_q - np.max(log_q))
            self.q[pair] /= np.sum(self.q[pair])
        self.qs.append(copy(self.q))

    def m_step(self):
        self.theta = np.zeros_like(self.theta, dtype=np.float)
        self.pkb = np.zeros_like(self.pkb, dtype=np.float)
        self.nu = np.array(self.nmut, dtype=np.float)

        for pair in self.mut:
            read, gene = pair
            kmers, inds = self.get_valid_kmers(pair)
            self.theta[inds] += self.q[pair]
            pkb_ind = self.get_pkb_ind(inds, read[2])
            self.pkb[inds, pkb_ind] += self.q[pair]
            self.nu[inds] += self.q[pair]

        self.theta = self.theta / (self.theta + self.nmut)
        self.theta[np.isnan(self.theta)] = 0
        self.thetas.append(copy(self.theta))

        self.pkb = (self.pkb.T / np.sum(self.pkb, axis=1)).T
        self.pkb[np.isnan(self.pkb)] = 0
        self.pkbs.append(copy(self.pkb))

        self.nu /= np.sum(self.nu)
        self.nu[np.isnan(self.nu)] = 0
        self.nus.append(copy(self.nu))


    def initial_params(self):
        self.nu = np.random.dirichlet([1] * self.all_k, size=1)[0]

        self.theta = np.sum(self.init_matr[:, 1:], axis=1) / np.sum(self.init_matr, axis=1)
        self.theta[np.isnan(self.theta)] = 0.05
        self.theta[self.theta < 1e-3] = 1e-3
        self.theta[self.theta > 1 - 1e-3] = 1 - 1e-3

        self.pkb = (self.init_matr[:, 1:].T / np.sum(self.init_matr[:, 1:], axis=1)).T
        self.pkb[np.isnan(self.pkb)] = 1
        self.pkb[self.pkb < 1e-3] = 1e-3
        self.pkb[self.pkb > 1 - 1e-3] = 1 - 1e-3
        self.pkb = (self.pkb.T / np.sum(self.pkb, axis=1)).T

        self.nus, self.thetas, self.pkbs = [copy(self.nu)], [copy(self.theta)], [copy(self.pkb)]
        self.qs = []
        self.e_step()
        self.m_step()

    def run(self):
        L_old, L = -np.inf, -np.inf
        iteration = 1
        self.Ls = []
        while (iteration < 100):
            iteration += 1
            self.e_step()
            L_old = self.get_L()
            self.m_step()
            L = self.get_L()
            print(L, L_old)
            assert L_old <= L
            self.Ls.append(L)
        # print(sorted(self.theta))
        # np.set_printoptions(threshold=np.nan)
        print(self.init_matr)
        print(np.round(self.theta, 3))
        print(np.round(self.thetas[0], 3))
        print(np.round(self.pkb, 3))
        print(np.round(self.pkbs[0], 3))

    def get_valid_kmers(self, pair):
        read, gene = pair
        kmers = [read[0], gene[0]]
        for i in xrange(1, len(read)):
            new_kmers = []
            for kmer in kmers:
                if i != 2:
                    new_kmers.append(kmer + read[i])
                new_kmers.append(kmer + gene[i])
            kmers = new_kmers
        kmers.sort()
        indices = map(lambda x: kmer_utils.kmer_index(x), kmers)
        return kmers, np.array(indices)

    def get_pkb_ind(self, I, B):
        B_ind = list("ACGT").index(B)
        centr_ind = (I % 4**3) / 4**2
        res = np.array([B_ind] * len(I))
        assert all(res != centr_ind)
        res[res > centr_ind] -= 1
        return res


    def get_L(self):
        assert np.sum(np.isnan(self.theta)) == 0
        assert np.sum(np.isnan(self.nu)) == 0
        assert np.sum(np.isnan(self.pkb)) == 0
        result = self.nmut.dot(np.nan_to_num(np.log(self.nu) + np.log(1 - self.theta)))
        for pair in self.mut:
            assert sum(np.isnan(self.q[pair])) == 0
            read, gene = pair
            kmers, inds = self.get_valid_kmers(pair)
            local_pkb = self.pkb[inds, self.get_pkb_ind(inds, read[2])]
            result += self.q[pair].dot(np.log(self.theta[inds]) + np.log(self.nu[inds]) + \
                                       np.log(local_pkb) - \
                                       np.log(self.q[pair]))
        return result


    def __init__(self, input_file):
        self.k = 5
        self.all_k = 4**self.k
        pinp = self.read_pinp(input_file)
        self.read_input(pinp)
        self.initial_params()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        type=str)
    params = parser.parse_args()
    em = EM(params.input)
    em.run()


if __name__ == "__main__":
    main()
