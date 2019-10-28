"""
Methods to perform decomposition with dynamic programming
"""

import itertools
from multitree import Tree

FOUTSIDE=0
FROOT=1
FINTERNAL=2
FSINGLE=3
FOREST='ORIS'
INFTY=1e100


def decompose(gt, st, weights=(1, 1000)):
    """

    :param gt: str or multitree.Tree object
    :param st: str or multitree.Tree object
    :param weights: tuple
    :return:
    """
    if not isinstance(gt, Tree):
        try:
            gt = Tree(gt)
        except:
            raise ValueError("Improper gene tree")
    if not isinstance(st, Tree):
        try:
            st = Tree(st)
        except:
            raise ValueError("Improper species tree")

    LOSS, ROOT = weights
    deltav = {}
    deltaupv = {}
    deltatreev = {}

    deltares = {}
    deltatreeres = {}
    deltaupres = {}

    def delta(g, s, inside, full):

        if (g, s) in deltav: return deltav[g, s]

        if g.leaf():
            if g.map == s:
                res = 0
            else:
                res = INFTY
            if full:
                deltares[g, s, inside] = []
        else:
            g1, g2 = g.c

            losscount = int(s.degree > 2) * LOSS * inside

            # speciation
            if not s.leaf():
                alpha = losscount + min(deltaup(g1, s1, 1, full) + deltaup(g2, s2, 1, full)
                                        for s1, s2 in itertools.permutations(s.c, 2))
            # print "HERE",g,s,alpha
            else:
                alpha = INFTY

            # hgt
            gamma = ROOT + min(
                deltatree(g1, g1.map, 0, full) + deltaup(g2, s, inside, full),
                deltatree(g2, g2.map, 0, full) + deltaup(g1, s, inside, full))
            res = min(alpha, gamma)
            # if s==st.root.c[1] and g==gt.root.c[0]:
            # 	print "HERE-G",g,s,gamma
            # 	print g1,g2,deltatree(g1,s,full),deltaup(g2,s,full)
            # 	print g1,g2,deltatree(g2,s,full),deltaup(g1,s,full)

            if full:
                deltares[g, s, inside] = []
                if not s.leaf() and alpha == res:
                    for s1, s2 in itertools.permutations(s.c, 2):
                        if res == losscount + deltaup(g1, s1, 1, full) + deltaup(g2, s2, 1, full):
                            for x, y in itertools.product(deltaupres[g1, s1, 1], deltaupres[g2, s2, 1]):
                                deltares[g, s, inside].append((x, y, "SPEC"))
                if gamma == res:
                    if res == ROOT + deltatree(g1, g1.map, 0, full) + deltaup(g2, s, inside, full):
                        for x, y in itertools.product(deltatreeres[g1, g1.map, 0], deltaupres[g2, s, inside]):
                            deltares[g, s, inside].append((x, y, "HGT1"))

                    if res == ROOT + deltatree(g2, g2.map, 0, full) + deltaup(g1, s, inside, full):
                        for y, x in itertools.product(deltatreeres[g2, g2.map, 0], deltaupres[g1, s, inside]):
                            deltares[g, s, inside].append((x, y, "HGT2"))

        deltav[g, s, inside] = res
        return res

    def deltaup(g, s, inside, full):
        if (g, s, inside) not in deltaupv:
            losscount = int(s.degree > 1) * LOSS * inside
            if s.leaf():
                res = delta(g, s, inside, full)
            else:
                res = min(delta(g, s, inside, full), losscount + min(deltaup(g, x, inside, full) for x in s.c))
            deltaupv[g, s, inside] = res
            if full:
                deltaupres[g, s, inside] = []
                if res == delta(g, s, inside, full): deltaupres[g, s, inside].append(s)
                if not s.leaf():
                    for x in s.c:
                        if res == losscount + deltaup(g, x, inside, full):
                            deltaupres[g, s, inside].extend(deltaupres[g, x, inside])
        return deltaupv[g, s, inside]

    def deltatree(g, s, inside, full):
        if (g, s, inside) not in deltatreev:
            deltatreev[g, s, inside] = res1 = min(delta(g, c, inside, full) for c in s.nodes())

            if full:
                deltatreeres[g, s, inside] = [c for c in s.nodes() if res1 == delta(g, c, inside, full)]

        return deltatreev[g, s, inside]

    def _gensolution(g, s, inside):
        if g.leaf():
            yield "%s" % g
        else:
            g1, g2 = g.c
            for x, y, tp in deltares[g, s, inside]:
                if tp == 'SPEC':
                    for t1, t2 in itertools.product(_gensolution(g1, x, 1), _gensolution(g2, y, 1)):
                        yield "(%s,%s) tp=%s" % (t1, t2, tp)
                elif tp == "HGT1":
                    for t1, t2 in itertools.product(_gensolution(g1, x, 0), _gensolution(g2, y, inside)):
                        yield "(%s linecolor='red',%s ) tp=%s" % (t1, t2, tp)
                elif tp == "HGT2":
                    for t1, t2 in itertools.product(_gensolution(g1, x, inside), _gensolution(g2, y, 0)):
                        yield "(%s,%s linecolor='green') tp=%s" % (t1, t2, tp)

    def _genforest(g, s, inside):
        if g.leaf():
            yield "%s" % g
        else:
            g1, g2 = g.c
            for x, y, tp in deltares[g, s, inside]:
                if tp == "SPEC":
                    for t1, t2 in itertools.product(_genforest(g1, x, 1), _genforest(g2, y, 1)):
                        yield "(%s,%s)" % (t1, t2)
                elif tp == "HGT1":
                    for t1, t2 in itertools.product(_genforest(g1, x, 0), _genforest(g2, y, inside)):
                        yield "(%s,%s,T)" % (t1, t2)
                else:
                    for t1, t2 in itertools.product(_genforest(g1, x, inside), _genforest(g2, y, 0)):
                        yield "(%s,%s,T)" % (t2, t1)

    def genforestroots(g, s):
        """
        Generates optimal sets of roots
        :param g: gene tree
        :param s: species tree
        :return:
        """
        def _genforestroots(g, s, inside):
            if g.leaf():
                yield g, []
            else:
                g1, g2 = g.c
                for x, y, tp in deltares[g, s, inside]:
                    if tp == "SPEC":
                        for (r1, t1), (r2, t2) in itertools.product(_genforestroots(g1, x, 1),
                                                                    _genforestroots(g2, y, 1)):
                            yield g, t1 + t2
                    elif tp == "HGT1":
                        for (r1, t1), (r2, t2) in itertools.product(_genforestroots(g1, x, 0),
                                                                    _genforestroots(g2, y, inside)):
                            yield r2, [r1] + t1 + t2
                    else:
                        for (r1, t1), (r2, t2) in itertools.product(_genforestroots(g1, x, inside),
                                                                    _genforestroots(g2, y, 0)):
                            yield r1, [r2] + t1 + t2

        for r, t in _genforestroots(g, s, 0): yield [r] + t

    gt.setlcamapping(st)

    for s in st.nodes:
        if s.leaf():
            s.degree = 0
        else:
            s.degree = len(s.c)

    deltatree(gt.root, st.root, 0, 1) + ROOT  # +the main tree

    # nodeid - postfix
    for s in deltatreeres[gt.root, st.root, 0]:
        for t in genforestroots(gt.root, s):
            yield tuple(i.nodeid for i in t)


