class MinHeap:
    import heapq

    def __init__(self, initial=None):
        if initial is None:
            initial = []
        self.heap = initial[:]
        self.heapq.heapify(self.heap)

    def push(self, item):
        self.heapq.heappush(self.heap, item)

    def pop(self):
        return self.heapq.heappop(self.heap)

    def pushpop(self, item):
        return self.heapq.heappushpop(self.heap, item)

    def pick(self):
        return self.heap[0]

    def __len__(self):
        return len(self.heap)


class MaxHeap:
    def __init__(self, initial=None):
        if initial is None:
            initial = []
        initial = [-item for item in initial]
        self.min_heap = MinHeap(initial)

    def push(self, item):
        self.min_heap.push(-item)

    def pop(self):
        return -self.min_heap.pop()

    def pushpop(self, item):
        return -self.min_heap.pushpop(-item)

    def pick(self):
        return -self.min_heap.pick()

    def __len__(self):
        return len(self.min_heap)


assert MaxHeap([1, 2, 3, 4]).pick() == 4
assert MinHeap([1, 2, 3, 4, 0]).pick() == 0


def cumulative_median(v, mean=lambda x, y: (x + y) / 2.):
    """O(n log n) two heaps algorithm"""

    max_heap = MaxHeap([float("-inf")])
    min_heap = MinHeap([float("+inf")])

    # Invariant: ALL max_heap <= ALL min_heap
    assert max_heap.pick() <= min_heap.pick()

    result = []
    for item in v:
        item = float(item)
        if item < max_heap.pick():
            item = max_heap.pushpop(item)
        elif item > min_heap.pick():
            item = min_heap.pushpop(item)

        assert max_heap.pick() <= item <= min_heap.pick()
        if len(max_heap) < len(min_heap):
            max_heap.push(item)
        else:
            min_heap.push(item)

        assert abs(len(max_heap) - len(min_heap)) <= 1
        if len(max_heap) > len(min_heap):
            median = max_heap.pick()
        elif len(max_heap) < len(min_heap):
            median = min_heap.pick()
        else:
            median = mean(max_heap.pick(), min_heap.pick())

        result.append(median)

    return result


def cumulative_mean(v):
    result = []
    S = 0
    for i, item in enumerate(v):
        S += item
        result.append(S / float(i + 1))
    return result


def cumulative_median_naive(v):
    import numpy as np
    v = np.array(v)
    result = []
    for i in xrange(1, len(v) + 1):
        result.append(np.median(v[:i]))
    return result


def reversed_cumulative_median(v):
    return cumulative_median(v[::-1])[::-1]


def reversed_cumulative_mean_numpy(v):
    import numpy as np

    return (np.cumsum(v[::-1]) / np.arange(1, len(v) + 1))[::-1]


def reversed_cumulative_mean(v):
    return cumulative_mean(v[::-1])[::-1]


def unique(x, y=None):
    if y is None:
        y = x

    assert len(x) == len(y)
    if len(x) == 0:
        return []

    result = [y[0]]
    for i in xrange(1, len(x)):
        if x[i] != x[i - 1]:
            result.append(y[i])

    return result

def how_many_greater_or_equal(x, limit):
    import bisect
    return len(x) - bisect.bisect_left(x, limit)

assert how_many_greater_or_equal([4, 4, 4, 4], 4) == 4
assert how_many_greater_or_equal([4], 4) == 1
assert how_many_greater_or_equal([], 4) == 0
assert how_many_greater_or_equal([1, 2, 3, 4, 5], 4) == 2
