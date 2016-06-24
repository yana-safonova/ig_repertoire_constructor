class KnuthIndex:

    def __init__(self, reads, tau=1, priority=None):
        from collections import defaultdict

        l = len(reads[0])

        for read in reads:
            assert(len(read) == l)

        piece_len = int(l / (tau + 1))
        piece_len_last = l - piece_len * tau

        piece_lens = [piece_len] * tau + [piece_len_last]

        if priority is None:
            priority = [0] * len(reads)
        self.priority = priority

        self.piece_lens = piece_lens
        self.piece_len = piece_len

        self.reads = reads
        self.sets_interval = []
        self.l = l
        self.tau = tau

        for j in range(self.tau + 1):
            sets = defaultdict(list)
            for i in range(len(reads)):
                substr = reads[i][j*piece_len:j*piece_len + piece_lens[j]]
                sets[substr].append(i)
            self.sets_interval.append(sets)

    def get_nn(self, q, nonn="None"):
        from Levenshtein import hamming

        assert(len(q) == self.l)

        set_for_check = set()

        for j in range(self.tau + 1):
            substr = q[j*self.piece_len:j*self.piece_len + self.piece_lens[j]]
            if substr in self.sets_interval[j]:
                set_for_check.update(self.sets_interval[j][substr])

        if len(set_for_check) == 0:
            return None

        i = min(set_for_check,
                key=lambda i: (hamming(self.reads[i], q), -self.priority[i]))

        result = self.reads[i]

        if hamming(result, q) > self.tau and nonn == "None":
            result = None

        return result


if __name__ == "__main__":
    reads = ["AAAA", "ACAA", "AAAC", "CAAA"]
    priority = [0, 1, 2, 3]
    priority = [0, -1, -2, -3]
    index = KnuthIndex(reads, 1, priority)

    print index.get_nn("CAAC")
    print index.get_nn("ACCC", nonn="Fatalistic")
    print index.get_nn("ACCC", nonn="None")
    print index.get_nn("CCCC", nonn="Fatalistic")
