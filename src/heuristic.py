"""
Methods to perform heuristic decomposition
"""

from generic import Decomposition


def assign_ranks(tree):
    """
    Assings taxonomic ranks (depths of a node) to a tree.
    Modifies the tree in situ.
    :param tree: ete2.Tree
    :return: None
    """
    T = tree
    for n in T.traverse(strategy="postorder"):
        if n.is_leaf():
            n.rank = 0
        else:
            n.rank = max(c.rank for c in n.children) + 1


def compute_mappings(G, S):
    """
    Computes both I and P mappings using Pawel's algorithm. Modifies G in situ.
    Leaf names in S are assumed to be unique.
    Leaf names in G are assumed to correspond to leaf names in S (no identifiers after underscore!)
    :param G: ete2.Tree
    :param S: ete2.Tree
    :return: None
    """
    try:
        d = S.rank
    except AttributeError:
        print("Species tree not initialized; Assign ranks before computing mappings.")
        raise

    # I mapping
    for g in G.traverse(strategy="postorder"):
        if g.is_leaf():
            g_species = g.name
            g_species = S.get_leaves_by_name(name=g_species)[0]
            g.M = g_species
            g.I = 0
            g.P = 0
            g.smap = g_species
        else:
            g_species = g.children[0].M.get_common_ancestor(g.children[1].M)
            g.M = g_species
            g.I = g_species.rank
            g.P = None

    # P mapping
    for s in S.traverse():
        s.lastvisited = None

    for i in range(0, d + 1):
        for g in G:
            if g.smap.rank == i:
                if g.smap.lastvisited is not None:
                    p = g.get_common_ancestor(g.smap.lastvisited)
                    if p.P is None:
                        p.P = i
                g.smap.lastvisited = g
                g.smap = g.smap.up


def unfold_tree(S):
    """
    Unfolds a species tree into evolutionary lineages.
    Returns a dictionary labelled by the leafs of S, containing the corresponding lineages
    as lists of tree nodes.
    The species tree needs to be labelled with taxonomic ranks by the assign_ranks() function.
    :param S: ete2.Tree
    :return: None
    """
    _, d = S.get_farthest_leaf()
    lineages = dict((s, []) for s in S)
    for s in S:
        p = s
        lge = lineages[s]
        while not p.is_root():
            pp = p.up
            lge.extend([p]*(pp.rank-p.rank))
            p = pp
        lge.append(p)
    assert all(len(lineages[l])==d+1 for l in lineages), "Error during tree unfolding. Please report this to the author"
    return lineages


def minimal_nodes(G):
    """
    Identifies the minimal nodes in G such that node.I != node.P or node.I == 1.
    Does not modify G.
    :param G: ete2.Tree object
        Gene tree decorated with I and P mapping values.
    :return: list
        A list of nodes from tree G.
        The elements of the list are references to nodes in G.
    """
    is_minimal = dict((g, False) for g in G.traverse())
    minimal_under = dict((g, False) for g in G.traverse())
    for g in G.traverse(strategy='postorder'):
        if g.is_leaf():
            continue
        minimal_under[g] = any(is_minimal[c] or minimal_under[c] for c in g.children)
        if not minimal_under[g]:
            if g.I == 0 or g.I > g.P:
                is_minimal[g] = True
    return [g for g in is_minimal if is_minimal[g]]


def cut_tree(g, lineages):
    """
    Finds an optimal cut under a node g with respect to the reference species tree,
    represented as a dictionary of evolutionary lineages.
    An optimal cut is chosen from edges adjacent to the children of g.
    If the node is embeddable, it is returned and a warning is issued.
    The function estimates the loss cost as shown in the paper.
    :param g: ete2.Tree object
        A node under which to cut (with I and P mappings).
    :param S: ete2.Tree object
        Reference species tree (with ranks assigned)
    :return: ete2.Tree
        The root of the new subtree (under the cut edge).
    """

    def check_division(leafset1, leafset2, r):
        """
        Checks if rank r divides evolutionary lineages of two sets of
        leaves from the gene tree
        (i.e. if the lineages of their species have disjoint taxa in rank r-1).
        :param r: int
            The rank, an integer from 0 to height of S.
        :param leafset: set
            A set of leaf nodes (ete3.Tree objects).
        :return: bool
        """
        if r == 0:
            # This would mean that we check an internal node with I == 0 after pruning,
            # so there is no embeddability.
            return False
        ancset1 = set([lineages[l.M][r - 1] for l in leafset1])
        ancset2 = set([lineages[l.M][r - 1] for l in leafset2])
        return ancset1.isdisjoint(ancset2)

    if g.I == g.P and g.I > 0:
        raise ValueError("The node is embeddable - nothing to cut!")
    elif g.is_leaf():
        raise ValueError("The node is a leaf - nothing to cut!")
    # Assigning candidate nodes for the root of the detached subtree
    dplc = g.children  # duplication-like candidates
    hgtc0 = dplc[0].children  # hgt-like candidates; either 0 or 2
    hgtc1 = dplc[1].children
    lc0 = len(hgtc0)
    lc1 = len(hgtc1)  # number of hgt-like candidates
    # dpl_emb = [True, True] #  g's children always satisfy embeddability condition
    dpl_emb = [dplc[0].I == g.I, dplc[1].I == g.I]  # alternative condition
    hgt0_emb = [False] * lc0
    hgt1_emb = [False] * lc1
    dpl_loss = [len(dplc[0].M) + len(dplc[1].M)] * 2
    hgt0_loss = [-1] * lc0
    hgt1_loss = [-1] * lc1
    dpl_size = [len(dplc[0]), len(dplc[1])]  # sizes of cut trees
    hgt0_size = [-1]*lc0
    hgt1_size = [-1]*lc1

    # checking embeddability and loss cost
    if hgtc0:
        leafset1 = set(dplc[1])
        # checking first element of hgtc0
        newM = dplc[1].M.get_common_ancestor(hgtc0[1].M)
        newI = newM.rank
        hgt0_emb[0] = check_division(set(hgtc0[1]), leafset1, newI)
        hgt0_loss[0] = len(hgtc0[0].M) + len(newM)
        hgt0_size[0] = len(hgtc0[0])
        # checking second element of hgtc0
        newM = dplc[1].M.get_common_ancestor(hgtc0[0].M)
        newI = newM.rank
        hgt0_emb[1] = check_division(set(hgtc0[0]), leafset1, newI)
        hgt0_loss[1] = len(hgtc0[1].M) + len(newM)
        hgt0_size[1] = len(hgtc0[1])

    if hgtc1:
        leafset1 = set(dplc[0])

        newM = dplc[0].M.get_common_ancestor(hgtc1[1].M)
        newI = newM.rank
        hgt1_emb[0] = check_division(set(hgtc1[1]), leafset1, newI)
        hgt1_loss[0] = len(hgtc1[0].M) + len(newM)
        hgt1_size[0] = len(hgtc1[0])

        newM = dplc[0].M.get_common_ancestor(hgtc1[0].M)
        newI = newM.rank
        hgt1_emb[1] = check_division(set(hgtc1[0]), leafset1, newI)
        hgt1_loss[1] = len(hgtc1[1].M) + len(newM)
        hgt1_size[1] = len(hgtc1[1])

    candidates = dplc + hgtc0 + hgtc1
    embeds = dpl_emb + hgt0_emb + hgt1_emb
    costs = dpl_loss + hgt0_loss + hgt1_loss
    sizes = dpl_size + hgt0_size + hgt1_size
    assert all(c >= 0 for c in costs), "Negative cost detected. Please report this to the authors"

    curr_candidate = candidates[0]
    curr_cost = costs[0]
    curr_size = sizes[0]
    # identify first embeddable candidate; used for the alternative embeddability condition
    # for the g's children
    for i, e in enumerate(embeds):
        if e:
            curr_candidate = candidates[i]
            curr_cost = costs[i]
            curr_size = sizes[i]
            break
    # identify the best candidate
    for i, c in enumerate(candidates[1:]):
        if embeds[i+1] and costs[i+1] < curr_cost:
            curr_candidate = c
            curr_cost = costs[i+1]
            curr_size = sizes[i+1]
        # elif embeds[i+1] and costs[i+1] == curr_cost and sizes[i+1] > curr_size:
        #     curr_candidate = c
        #     curr_cost = costs[i + 1]
        #     curr_size = sizes[i + 1]

    return curr_candidate

    #
    # for i in range(lc1):
    #     # if hgtc1 not empty, then it has two elements; 1-i is the other one
    #     leafset1 = set(hgtc1[1 - i])
    #     newM = dplc[1].M.get_common_ancestor(hgtc1[1 - i].M)  # updated M mapping after removal of hgtc1[i]
    #     newI = newM.rank
    #     embeddable[2 + i] = check_division(leafset1, leafset2, newI)
    #     loss_cost[2 + i] = len(hgtc1[i].M) + len(newM)  # turned out that investigating size of G is unneccessary
    #     new_tree_size[2 + i] = len(hgtc1[i])
    #
    # leafset2 = set(dplc[0])
    # for i in range(lc2):
    #     leafset1 = set(hgtc2[1 - i])
    #     newM = dplc[0].M.get_common_ancestor(hgtc2[1 - i].M)
    #     newI = newM.rank
    #     embeddable[2 + lc1 + i] = check_division(leafset1, leafset2, newI)
    #     loss_cost[2 + lc1 + i] = len(hgtc2[i].M) + len(newM)
    #     new_tree_size[2 + lc1 + i] = len(hgtc2[i])

    # # not very efficient, but oh well
    # min_cost = min([loss_cost[i] for i in range(2 + lc1 + lc2) if embeddable[i]])
    # # embeddable & min cost:
    # good_cut = [loss_cost[i] == min_cost and embeddable[i] for i in range(2 + lc1 + lc2)]
    # target_size = min([new_tree_size[i] for i in range(2 + lc1 + lc2) if good_cut[i]])
    # # embeddable, min cost & min tree:
    # cut_id = [i for i in range(2 + lc1 + lc2) if good_cut[i] and new_tree_size[i] == target_size][0]
    #
    # if cut_id < 2:
    #     newroot = dplc[cut_id]
    # elif cut_id < 2 + lc1:
    #     newroot = hgtc1[cut_id - 2]
    # elif cut_id < 2 + lc1 + lc2:
    #     newroot = hgtc2[cut_id - lc1 - 2]
    # else:
    #     raise RuntimeError("Wrong rootnode id!")
    # return newroot


def decompose(gene_tree, species_tree):
    """
    Performs the decomposition of the gene tree with respect to the species tree.
    Returns the forest of subtrees as a list.
    :param gene_tree: ete2.Tree object
    :param species_tree:  ete2.Tree object
        The reference species tree. The leaf names need to be unique. The tree needs to be
        ranked, i.e. each node needs to have an integer 'rank' attribute. Ranks
        can be assigned e.g. by running the function assign_ranks(species_tree) from this package.
    :return: Decomposition object
        A Decomposition of the gene tree.
    """
    G = gene_tree.copy()
    S = species_tree.copy()  # the tree will be modified in assign_ranks function

    # checking ranks of S
    for s in S.traverse():
        if not hasattr(s, 'rank'):
            raise ValueError("""Species tree is not ranked.\nPlease provide ranks e.g. by running assign_ranks(species_tree).""")
        elif s.is_leaf() and s.rank != 0:
            raise ValueError("Leaf of a species tree not ranked as 0")

##    # validate the labelling of leafs of S
##    snames = set()
##    for s in S:
##        if s.name in snames:
##            raise ValueError("Node names of the species tree are not unique!")
##        else:
##            snames.add(s.name)

##    # names of leaves of G are stripped from the identifiers; original names are not used later,
##    # because original gene tree will be returned
##    for g in G:
##        newname = g.name.split('_')[0]
##        if newname not in snames:
##            raise ValueError("The gene %s does not correspond to any species!" % newname)
##        g.name = newname

    for i, g in enumerate(G.traverse(strategy="postorder")):
        g.nid = i

    # initialize data
    lineages = unfold_tree(S)

    # decompose
    improper = [G]  # list of minimal non-embeddable nodes
    roots = []  # decomposition forest
    while improper:
        compute_mappings(G, S)
        improper = minimal_nodes(G)
        subtrees = [cut_tree(g, lineages) for g in improper]
        roots += [s.nid for s in subtrees]
        for s in subtrees:
            s.detach()
        # Pruning G. This is time costly, and ideally should be get rid of,
        # but is required in the current implementation.
        while len(G.children) == 1:
            G = G.children[0]  # moving down the root, because pruning preserves it
        G.prune(G)  # removing nodes with single children
    roots.append(G.nid)
    return roots


if __name__ == "__main__":
    from time import time

    #G = ete3.Tree('((((((((S21_1:1,(S22_3:1,(S22_9:1,S22_10:1)1:1)1:1)1:1,((S24_1:1,(S25_1:1,S26_1:1)1:1)1:1,((S29_9:1,(S35_14:1,S6_4:1)1:1)1:1,(S0_5:1,S0_6:1)1:1)1:1)1:1)1:1,S32_1:1)1:1,(((((S21_7:1,(S21_11:1,(S22_15:1,S38_7:1)1:1)1:1)1:1,S21_4:1)1:1,(S22_5:1,S22_6:1)1:1)1:1,(((S25_3:1,S26_3:1)1:1,S29_6:1)1:1,(S25_2:1,(S49_3:1,S26_4:1)1:1)1:1)1:1)1:1,(S32_3:1,S32_4:1)1:1)1:1)1:1,((S34_4:1,S35_2:1)1:1,S35_1:1)1:1)1:1,(((S38_3:1,((S63_3:1,(S64_3:1,S65_3:1)1:1)1:1,(S24_3:1,S43_4:1)1:1)1:1)1:1,((S39_2:1,S39_3:1)1:1,((S40_1:1,S41_1:1)1:1,(S43_3:1,S44_1:1)1:1)1:1)1:1)1:1,(S49_1:1,((S50_2:1,S50_3:1)1:1,(S51_2:1,S29_2:1)1:1)1:1)1:1)1:1)1:1,(((S56_1:1,S57_1:1)1:1,((S60_2:1,S60_3:1)1:1,S34_2:1)1:1)1:1,((S56_2:1,(S29_1:1,S57_3:1)1:1)1:1,((S63_1:1,(S64_1:1,S65_1:1)1:1)1:1,(S63_2:1,(S64_2:1,S65_2:1)1:1)1:1)1:1)1:1)1:1)1:1,(S19_2:1,(S19_3:1,(((((((S21_9:1,S21_10:1)1:1,(S35_10:1,(((S25_5:1,S26_9:1)1:1,S22_17:1)1:1,(S22_18:1,S22_19:1)1:1)1:1)1:1)1:1,(((S43_6:1,S24_5:1)1:1,(S25_4:1,(S26_7:1,S26_8:1)1:1)1:1)1:1,S29_7:1)1:1)1:1,(S32_6:1,S32_7:1)1:1)1:1,(((S21_6:1,(S2_1:1,S22_12:1)1:1)1:1,(S29_8:1,(S56_4:1,(S57_8:1,S57_9:1)1:1)1:1)1:1)1:1,((S35_11:1,S35_12:1)1:1,(S61_1:1,S35_13:1)1:1)1:1)1:1)1:1,((S49_2:1,(S60_4:1,(S50_4:1,S51_3:1)1:1)1:1)1:1,(S39_4:1,((S40_2:1,(S41_4:1,S26_6:1)1:1)1:1,((S21_12:1,S43_7:1)1:1,S44_2:1)1:1)1:1)1:1)1:1)1:1,(((((S35_7:1,S35_8:1)1:1,S0_3:1)1:1,((S1_2:1,(S2_3:1,(S2_8:1,S2_9:1)1:1)1:1)1:1,S4_1:1)1:1)1:1,(((((S1_3:1,(S2_10:1,S1_6:1)1:1)1:1,(S34_6:1,S61_2:1)1:1)1:1,(((S0_7:1,S1_7:1)1:1,(S2_11:1,S60_5:1)1:1)1:1,S4_3:1)1:1)1:1,(((S1_5:1,(S39_5:1,(S2_15:1,S2_16:1)1:1)1:1)1:1,(S4_6:1,S4_7:1)1:1)1:1,S6_3:1)1:1)1:1,((S41_2:1,S8_1:1)1:1,(((S11_5:1,S11_6:1)1:1,(S12_2:1,S13_2:1)1:1)1:1,(S11_4:1,(((S12_4:1,S13_4:1)1:1,(S12_5:1,S13_5:1)1:1)1:1,(S12_3:1,S13_3:1)1:1)1:1)1:1)1:1)1:1)1:1)1:1,((((S11_2:1,S1_1:1)1:1,(S12_1:1,S13_1:1)1:1)1:1,(S56_3:1,(S57_5:1,S57_6:1)1:1)1:1)1:1,S19_4:1)1:1)1:1)1:1)1:1)1:1);')
    #S = ete3.Tree('(((S0:1,((((S1:1,S2:1)1:1,S4:1)1:1,S6:1)1:1,((S8:1,S9:1)1:1,(S11:1,(S12:1,S13:1)1:1)1:1)1:1)1:1)1:1,S19:1)1:1,(((((((S21:1,S22:1)1:1,((S24:1,(S25:1,S26:1)1:1)1:1,S29:1)1:1)1:1,S32:1)1:1,(S34:1,S35:1)1:1)1:1,((S38:1,(S39:1,((S40:1,S41:1)1:1,(S43:1,S44:1)1:1)1:1)1:1)1:1,(S49:1,(S50:1,S51:1)1:1)1:1)1:1)1:1,(S56:1,S57:1)1:1)1:1,(((S60:1,S61:1)1:1,(S63:1,(S64:1,S65:1)1:1)1:1)1:1,S69:1)1:1)1:1);')
    #assign_ranks(S)
    #roots = [1, 2, 14, 17, 24, 25, 27, 34, 40, 44, 50, 57, 61, 66, 67, 72, 85, 89, 98, 105, 113, 124, 125, 128, 131, 134, 140, 144, 152, 157, 163, 166, 168, 169, 172, 179, 188, 191, 200, 206, 207, 213, 217, 220, 221, 223, 224, 228, 235, 236, 241, 244, 248, 251, 257, 261, 274, 276, 281, 283, 288]

    G = ete3.Tree('(((a, b), ((a, c), (c, d))), ((b, d), f));')
    S = ete3.Tree('(a, b, c, d, e, f);')
    assign_ranks(S)
    #roots = [0, 1, 5, 8, 13, 14]

    # TD = Decomposition(G, S, roots)
    # TD.assign_topological_ranks()
    # TD.forest()
    # TD.locus_tree()
    #G = ete3.Tree('(a, (a,b));', format=8)
    #S = ete3.Tree('(a, b);', format=8)
    
    s = time()
    R = decompose(G, S)
    e = time()
    print("Decomposed in %.02f seconds" % (e - s))
    D = Decomposition(G, S, R)
    s = time()
    F = D.forest()
    e = time()
    print("Forest obtained in %.02f seconds" % (e - s))
    s = time()
    L = D.locus_tree()
    e = time()
    print("Locus tree obtained in %.02f seconds" % (e - s))
    D._assert_sources()
    print(L.get_ascii(attributes=['name', 'nid', 'source', 'type', 'I', 'P']))
    clstrs = [g for g in L.traverse() if g.inspect and (g.is_root() or not g.up.inspect)]

    # for i in range(1000):
    #     if not i % 100: print i
    #     S = ete2.Tree()
    #     S.populate(20)
    #     G = ete2.Tree()
    #     G.populate(20, names_library=[s.name for s in S], reuse_names=True)
    #     D = decompose(G, S)
    #     F = D.forest()
    #     L = D.locus_tree()
    #     D._assert_sources()
        # lnms = [set([l.name for l in L.search_nodes(locus_id=i) if l.is_leaf()]) for i in set(l.locus_id for l in L)]
        # fnms = [set([l.name for l in f]) for f in F]

    #     assert all([l in fnms for l in lnms])
    #     assert all([l in lnms for l in fnms])
