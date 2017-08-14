from igquast_impl import *


def read_graph(fname):
    with smart_open(fname) as f:
        header = f.next()
        cons_len, E, FORMAT = map(int, header.strip().split())

        g = [[] for _ in xrange(cons_len)]

        for i, l in enumerate(f):
            l = l.strip().split()
            l = map(int, l)
            l = l[1:]
            g[i] = [n - 1 for n in l[0::2]]

        return g


def conn_components(g, top=10):
    N = len(g)
    comps = [None] * N
    def go(i, comp):
        if comps[i] is None:
            comps[i] = comp
            for j in g[i]:
                go(i, comp)

    for i in xrange(N):
        go(i, i)

    from collections import defaultdict
    dd = defaultdict(list)
    for i, comp in enumerate(comps):
        dd[comp].append(i)

    lst = dd.values()

    lst.sort(key=len, reverse=True)

    lst = lst[:top]
    return lst


def get_cleaned_nums(mapfile, nums):
    nums = set(nums)
    result = []
    with smart_open(mapfile) as f:
        for cleaned_num, line in enumerate(f):
            compressed_num = int(line.strip())
            if compressed_num in nums:
                result.append(cleaned_num)

    return result


def get_read_ids(readfile, nums):
    nums = set(nums)
    result = []
    with smart_open(readfile) as f:
        for i, record in enumerate(SeqIO.parse(f, idFormatByFileName(readfile))):
            if i in nums:
                result.append(str(record.description))
    return result


def subset_reads_by_id(input, ids, output, prob=0.05):
    import random
    random.seed(0)
    ids = set(ids)
    n = 0
    with smart_open(input) as fi, smart_open(output, "w") as fo:
        for record in SeqIO.parse(fi, idFormatByFileName(input)):
            if str(record.description) in ids or random.uniform(0., 1.) < prob:
                SeqIO.write(record, fo, idFormatByFileName(output))
                n += 1
    print "%d reads written into %s" % (n, output)


if __name__ == "__main__":
    dir = "../12/"
    graph_name = dir + "/sw.graph"
    map_name = dir + "/cleaned_compressed_map.txt"
    cleaned_reads = dir + "/vj_finder/cleaned_reads.fa"
    filtered_reads = dir + "/vj_finder/filtered_reads.fa"
    initial_reads = dir + "/merged_reads.fastq"

    g = read_graph(graph_name)
    conn_comps = conn_components(g, 10)
    nums = sum(conn_comps, [])
    print len(nums), "#compressed reads in top connectivity components"

    cleaned_nums = get_cleaned_nums(map_name, nums)
    print len(cleaned_nums), "#cleaned reads"
    ids = get_read_ids(cleaned_reads, cleaned_nums)
    print len(ids), "#ids"
    subset_reads_by_id(initial_reads, ids, "workshop_dataset.fq.gz")
