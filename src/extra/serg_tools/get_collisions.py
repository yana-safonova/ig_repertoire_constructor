import sys
from collections import defaultdict


def main():
    rcm = sys.argv[1]
    umi_to_cnt = defaultdict(int)
    for line in open(rcm):
        if line.startswith('original'):
            umi = line.split('\t')[0].split(':')[-1]
            umi_to_cnt[umi] += 1
    print len(umi_to_cnt)
    for key, value in umi_to_cnt.iteritems():
        if value > 1:
            print key


if __name__ == '__main__':
    main()
