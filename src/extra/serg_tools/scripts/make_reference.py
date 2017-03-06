import sys

# Converts simulated sequences to compressed
first = True
for line in sys.stdin:
    t = line.split("_")
    if len(t) == 1:
        if first:
            print t[0].strip()
    else:
        first = int(t[-1]) == 1
        if not first:
            continue
        print ">cluster___%d___size___%d" % (int(t[1]) - 1, int(t[3]))
