import os
import shutil
import sys


def main():
    log_path = sys.argv[1]
    if os.path.exists("bad"):
        shutil.rmtree("bad")
    os.makedirs("bad")
    with open(log_path, "rb") as log_file:
        for line in log_file:
            if "fasta" in line:
                line = line[:-1]
                shutil.copy2(os.path.join("all", line[:-len(".aln")]), "bad")
                shutil.copy2(os.path.join("clustal", line), "bad")


if __name__ == '__main__':
    main()
