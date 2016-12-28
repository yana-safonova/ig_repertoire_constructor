#!/usr/bin/env python2

from Bio import SeqIO
import os
import os.path
import re
import glob
import argparse
import sys
import sqlite3


def initDB(cursor):
    cursor.execute("create table reads ( \
                   id integer primary key not null, \
                   seqid string, \
                   read string, \
                   quality string, \
                   component integer, \
                   cluster integer, \
                   barcode string)")

from read_barcode_splitter import extract_barcode


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split source fastq file to hamming graph connectivity components")
    parser.add_argument("-s", "--reads",
                        type=str,
                        required=True,
                        help="Reads file in FASTQ format")
    parser.add_argument("-i", "--igrc-output-dir",
                        type=str,
                        required=True,
                        help="Ig_rep_constructor output directory")
    parser.add_argument("-o", "--output-dir",
                        type=str,
                        required=True,
                        help="Output directory")
    parser.add_argument("-f", "--force-remove-existing-output-dir",
                        action="store_true",
                        help="Write to existing output dir (default no)")
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help="Be verbose (default no)")
    parser.add_argument("-C", "--do-not-cut-reads",
                        action="store_true",
                        help="Don't cut reads, preserve original ones (default no)")
    parser.add_argument("-S", "--sort",
                        action="store_true",
                        help="Sort outputted reads by clique index (default no)")
    parser.add_argument("--sqlite", "-Q",
                        type=str,
                        help="output database file (in SQLite format)")
    parser.add_argument("--unclussified-output", "-B",
                        type=str,
                        help="output file for unclussified items")

    args = parser.parse_args()

    if args.sqlite is not None:
        if os.path.isfile(args.sqlite):
            os.unlink(args.sqlite)
        conn = sqlite3.connect(args.sqlite, isolation_level="EXCLUSIVE")
        initDB(conn)

    with open(args.reads, "rU") as fh:
        reads_hash_table = {record.id: record for record in SeqIO.parse(fh, "fastq")}

    if args.verbose:
        print("Loaded %d reads" % len(reads_hash_table))

    # Create output dir if needed
    try:
        os.mkdir(args.output_dir)
    except Exception as ex:
        if not args.force_remove_existing_output_dir:
            print(ex)
            sys.exit(1)
        else:
            print("Writing to existing output dir...")

    if not os.path.exists(args.igrc_output_dir):
        print("Passed Ig RC output dir don't exist")
        sys.exit(1)

    if not os.path.exists("%s/dense_subgraphs" % args.igrc_output_dir):
        print("Passed Ig RC output dir don't contain dense_subgraphs")
        print("Maybe you forget --output-dense-sgraphs option?")
        sys.exit(1)

    used_ids = set()

    for decomposition in glob.glob("%s/dense_subgraphs/*.txt" % args.igrc_output_dir):
        basename = os.path.basename(decomposition)

        m = re.match(r"dense_sgraph_decomposition_(\d+)_size_(\d+)\.txt",
                     basename)
        concomp_id = m.groups()[0] if m else None

        out_file_cliques = "%s/%s" % (args.output_dir, re.sub("txt$", "cliques.txt", basename))
        out_file_fastq = "%s/%s" % (args.output_dir, re.sub("txt$", "fastq", basename))
        with open(decomposition, "r") as dec_fh:
            res = []
            for line in dec_fh:
                cluster, seq_id, align = line.split()
                cluster, align = int(cluster), int(align)
                k = re.sub(r"_SUBSTR\(\d*,\d*\)", "", seq_id)

                used_ids.add(k)

                if (k != seq_id):
                    if args.verbose:
                        print("Substitute %s -> ..." % seq_id)
                    m = re.match(r".*_SUBSTR\((\d*),\d*\)", seq_id)
                    d = int(m.groups()[0])
                    align += d
                seq = reads_hash_table[k]
                if not args.do_not_cut_reads:
                    # Cut read head
                    seq = seq[align:]
                res.append((cluster, seq))

            min_len = min([len(seq_) for cluster_, seq_ in res])
            if not args.do_not_cut_reads:
                # Trim reads to minimal length
                res = [(cluster_, seq_[:min_len]) for cluster_, seq_ in res]

            if args.sort:
                # sort by cluster_id
                res.sort(key=lambda x: x[0])

            with open(out_file_cliques, "w") as out_fh:
                out_fh.write(" ".join([str(cluster_) for cluster_, record_ in res]))

            with open(out_file_fastq, "w") as out_fh:
                records = [record_ for cluster_, record_ in res]
                SeqIO.write(records, out_fh, "fastq")

            if args.sqlite is not None:
                for cluster, record in res:
                    conn.execute("insert into reads (seqid, read, quality, component, cluster, barcode) \
                                 values (?, ?, ?, ?, ?, ?)", \
                                 (str(record.id),
                                  str(record.seq),
                                  None,  # TODO extract quality
                                  concomp_id,
                                  cluster,
                                  extract_barcode(record.id)))
                conn.commit()

            if args.verbose:
                print("Processed file %s, %d lines, min = %d" % (decomposition,
                                                                 len(res),
                                                                 min_len))
    unclussified_reads = [record for id, record in reads_hash_table.iteritems() if id not in used_ids]

    if args.unclussified_output is not None:
        with open(args.unclussified_output, "w") as out_fh:
            SeqIO.write(unclussified_reads, out_fh, "fastq")

    if args.sqlite is not None:
        for record in unclussified_reads:
            conn.execute("insert into reads (seqid, read, quality, component, cluster, barcode) \
                            values (?, ?, ?, ?, ?, ?)", \
                            (str(record.id),
                            str(record.seq),
                            None,  # TODO extract quality
                            None,
                            None,
                            extract_barcode(record.id)))
        conn.commit()

    if args.sqlite is not None:
        conn.commit()

        conn.execute("create index if not exists read_seqid_index on reads(seqid)")
        conn.execute("create index if not exists read_component_index on reads(component)")
        conn.execute("create index if not exists read_component_cluster_index on reads(component, cluster)")
        conn.execute("create index if not exists read_barcode_index on reads(barcode)")
        conn.commit()

        print("TABLE `reads` indexed")

    if args.sqlite is not None:
        conn.commit()
        conn.close()
