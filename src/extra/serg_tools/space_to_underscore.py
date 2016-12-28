import sys


def main():
    inf = sys.argv[1]
    outf = sys.argv[2]
    with open(outf, "w") as outh:
        for line in open(inf):
            outh.write(line.replace(' ', '_'))

if __name__ == '__main__':
    main()
