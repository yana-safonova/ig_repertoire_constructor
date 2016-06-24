
# TODO factor in to the separate module
def count_matches(s, reference, align=0):
    if align < 0:
        s, reference = reference, s
        align = -align

    reference, s = str(reference), str(s)
    reference = reference[align:]

    return sum(c1 == c2 for c1, c2 in zip(s, reference))


