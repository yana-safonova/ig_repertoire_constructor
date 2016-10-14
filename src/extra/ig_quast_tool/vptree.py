from collections import deque
import random
import numpy as np


class BKTree:
    def __init__(self, points, dist_fn):
        from collections import defaultdict

        self.children = None
        self.dist_fn = dist_fn

        # choose a better vantage point selection process
        self.vp = points.pop(random.randrange(len(points)))

        if len(points) < 1:
            return

        # choose division boundary at median of distances
        distances = [self.dist_fn(self.vp, p) for p in points]

        subsamples = defaultdict(list)
        for d, p in zip(distances, points):
            subsamples[d].append(p)

        self.children = {d: BKTree(subsample, dist_fn)
                         for d, subsample in subsamples.iteritems()}

    def is_leaf(self):
        return self.children is None

    def get_all_in_range(self, q, tau):
        """
        find all points within a given radius of point q

        :param self: bk-tree
        :param q: a query point
        :param tau: the maximum distance from point q
        """

        # buffer for nearest neightbors
        neighbors = []

        # list of nodes ot visit
        visit_stack = deque([self])

        while len(visit_stack) > 0:
            node = visit_stack.popleft()
            if node is None:
                continue

            d = self.dist_fn(q, node.vp)
            if d <= tau:
                neighbors.append((d, node.vp))

            if node.is_leaf():
                continue

            for dd in range(d - tau, d + tau + 1):
                if dd in node.children:
                    visit_stack.append(node.children[dd])

        return neighbors


class VPTree(object):
    """
    An efficient data structure to perform nearest-neighbor
    search.
    """

    def __init__(self, points, dist_fn=None):
        self.left = None
        self.right = None
        self.mu = None
        self.dist_fn = dist_fn if dist_fn is not None else l2

        # choose a better vantage point selection process
        self.vp = points.pop(random.randrange(len(points)))

        if len(points) < 1:
            return

        # choose division boundary at median of distances
        distances = [self.dist_fn(self.vp, p) for p in points]
        self.mu = np.median(distances)

        left_points = []  # all points inside mu radius
        right_points = []  # all points outside mu radius
        for d, p in zip(distances, points):
            if d >= self.mu:
                right_points.append(p)
            else:
                left_points.append(p)

        del distances
        del d
        del p
        del points

        if len(left_points) > 0:
            self.left = VPTree(points=left_points, dist_fn=self.dist_fn)

        if len(right_points) > 0:
            self.right = VPTree(points=right_points, dist_fn=self.dist_fn)

    def is_leaf(self):
        return (self.left is None) and (self.right is None)

    def get_nearest_neighbors(self, q, k=1):
        """
        find k nearest neighbor(s) of q

        :param self:  vp-tree
        :param q: a query point
        :param k: number of nearest neighbors

        """

        # buffer for nearest neightbors
        neighbors = PriorityQueue(k)

        # list of nodes ot visit
        visit_stack = deque([self])

        # distance of n-nearest neighbors so far
        tau = np.inf

        while len(visit_stack) > 0:
            node = visit_stack.popleft()
            if node is None:
                continue

            d = self.dist_fn(q, node.vp)
            if d < tau:
                neighbors.push(d, node.vp)
                tau, _ = neighbors.queue[-1]

            if node.is_leaf():
                continue

            if d < node.mu:
                if d < node.mu + tau:
                    visit_stack.append(node.left)
                if d >= node.mu - tau:
                    visit_stack.append(node.right)
            else:
                if d >= node.mu - tau:
                    visit_stack.append(node.right)
                if d < node.mu + tau:
                    visit_stack.append(node.left)
        return neighbors.queue

    def get_all_in_range(self, q, tau):
        """
        find all points within a given radius of point q

        :param self: vp-tree
        :param q: a query point
        :param tau: the maximum distance from point q
        """

        # buffer for nearest neightbors
        neighbors = []

        # list of nodes ot visit
        visit_stack = deque([self])

        while len(visit_stack) > 0:
            node = visit_stack.popleft()
            if node is None:
                continue

            d = self.dist_fn(q, node.vp)
            if d <= tau:
                neighbors.append((d, node.vp))

            if node.is_leaf():
                continue

            if d + tau < node.mu:
                visit_stack.append(node.left)
            elif d - tau >= node.mu:
                visit_stack.append(node.right)
            else:
                visit_stack.append(node.left)
                visit_stack.append(node.right)

        return neighbors


class PriorityQueue(object):
    def __init__(self, size=None):
        self.queue = []
        self.size = size

    def push(self, priority, item):
        self.queue.append((priority, item))
        self.queue.sort()
        if self.size is not None and len(self.queue) > self.size:
            self.queue.pop()


def l2(p1, p2):
    return np.sqrt(np.sum(np.power(p2 - p1, 2)))




calls = 0
def lev_dist(p1, p2):
    from Levenshtein import distance

    global calls
    calls += 1

    return distance(p1, p2)

if __name__ == '__main__':
    import resource, sys
    resource.setrlimit(resource.RLIMIT_STACK, (2**29,-1))
    sys.setrecursionlimit(10**8)


    from Bio import SeqIO

    print("Data reading...")
    with open("constructed_repertoire.clusters.fa", "rU") as fh:
        reads = [str(read.seq) for read in SeqIO.parse(fh, "fasta")]

    reads = reads[:1000]

    print("Data loaded (N = %d)" % len(reads))

    print("VP-tree construction...")
    tree = VPTree(reads, lev_dist)
    print("VP-tree constructed")

    S = 0
    for i, q in enumerate(reads):
        # k_neighbors = tree.get_nearest_neighbors(q, 3)
        eps_neib = tree.get_all_in_range(q, 4)
        # S += len(eps_neib) + len(k_neighbors)
        print len(eps_neib)
        # if i % 100 == 0:
        #     print "Iteration %d" % i

        # print "query:"
        # print "\t", q
        # print "nearest neighbors: "
        # for d, n in neighbors:
        #     print "\t", n
    print S
    print calls
    print len(reads) ** 2 / 2
