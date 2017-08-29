import numpy as np
from Bio import SeqIO

def ones_runs(a):
    # Create an array that is 1 where a is 0, and pad each end with an extra 0.
    iszero = np.concatenate(([0], np.equal(a, 1).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where absdiff is 1.
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges


def count_distances(filename):
    records = list(SeqIO.parse(filename, "fasta"))
    genes = [str(record.seq) for record in records[1::2]]
    reads = [str(record.seq) for record in records[::2]]
    n_muts, distances, durations = [], [], []
    for gene, read in zip(genes, reads):
        if '-' in gene or '-' in read:
            continue
        pos = [i for i in xrange(len(gene)) if gene[i] != read[i]]
        pos = np.array(pos)
        diff = np.diff(pos)
        distances += list(diff)
        ones = ones_runs(diff)
        ones = ones[:,1] - ones[:,0]
        durations += list(ones)
        n_muts.append(len(pos))

    return distances, durations, n_muts


def main():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--input", "-i", required=True)

    params = parser.parse_args()
    with open(params.input, 'r') as f:
        input_files = [line.strip('\n') for line in f.readlines()]
    all_distances, all_durations, all_n_muts = [], [], []
    for filename in input_files:
        print(filename)
        distances, durations, n_muts = count_distances(filename)
        all_distances += distances
        all_durations += durations
        all_n_muts += n_muts
    print("\nDistances")
    from collections import Counter
    cnt = Counter(all_distances)
    print(cnt)
    print("\nDurations")
    dur_1 = len(filter(lambda x: x == 1, all_durations))
    print("Single mutations : %d/%d" % (dur_1, len(all_durations)))
    cnt = Counter(all_durations)
    print(cnt)
    print("\nN_MUTS")
    cnt = Counter(all_n_muts)
    print(cnt)


if __name__ == "__main__":
    main()
