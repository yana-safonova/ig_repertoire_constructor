import sys

size = 0
for line in sys.stdin:
    t = line.split("___")
    if len(t) > 1:
        size += int(t[-1])
print size
