import ete2
from warnings import warn


class Decomposition(object):
    def __init__(self, gene_tree, species_tree, roots):
        """
        A locus decomposition of a gene tree with respect to a given species tree.
        :param gene_tree: ete2.Tree
        :param species_tree: ete2.Tree
        :param roots: list
            A list of identifiers (in post-order numbering) of the roots of trees from the decomposition forest.
        """
        self.G = gene_tree
        for i, g in enumerate(self.G.traverse(strategy="postorder")):
            g.nid = i
        self.G_labelled = False # whether G has been labelled with raw I/P mappings
        self.S = species_tree
        self.roots = sorted(roots)
        self.F = []  # forest
        self.L = None  # locus tree
        # self._assign_locus_ids()

    def _assign_locus_ids(self):
        """
        Assigns IDs of loci to nodes of the gene tree.
        :return:
        """
        locus = 0
        loci_to_process = [self.L]
        while loci_to_process:
            root = loci_to_process.pop()
            root.locus_id = locus
            for n in root.iter_descendants(strategy="postorder",
                                           is_leaf_fn=lambda x: x.source and x != root):
                # this will initially assign an improper locus id to a source node,
                # but it will be overriden in the next iteration:
                n.locus_id = locus
                if n.source:
                    loci_to_process.append(n)
            locus += 1

    def _get_node_types(self):
        """
        Classifies nodes of the locus tree (self.L) as either root, internal, or outside nodes.
        :return: dict
            A dictionary mapping gene tree's nodes' postorder IDs to their types.
        """
        nodetype = dict((g.nid, "internal") for g in self.L.traverse())
        for r in self.roots:
            nodetype[r] = "root"
        for g in self.L.traverse(strategy="postorder"):
            if nodetype[g.nid] == "root" and not g.is_root():
                nodetype[g.up.nid] = "outside"  # junction nodes
        for g in self.L.traverse(strategy="postorder"):
            if not g.is_leaf() and all(nodetype[c.nid] == "outside" for c in g.children):
                nodetype[g.nid] = "outside"  # nodes above junction nodes
        return nodetype

    def _compute_adjusted_mappings(self):
        """
        Computes the I/P mappings based on neighbouring loci.
        For any node, when computing it's mapping values, the tree is temporarily pruned
        to contain only the loci of the node's children.
        :return: ete2.Tree object
            A tree with I and P mapping values added as attributes to each node.
            The value of mapping X is stored as node.X attribute.
        """
        S = self.S
        L = self.L
        try:
            d = S.rank
        except AttributeError:
            print("Species tree not initialized; Assign ranks before computing mappings.")
            raise

        try:
            L.locus_id
        except AttributeError:
            print("Locus tree not initialized; Assign locus IDs before computing mappings.")
            raise

        # I mapping
        pureM = dict((g.nid, S) for g in L.traverse(strategy="postorder"))  # I from single loci
        combinedM = pureM.copy()  # I at junction nodes from neighbouring loci
        smap = dict((g, S) for g in L)
        for g in L.traverse(strategy="postorder"):
            if g.is_leaf():
                g_species = g.name.split('_')[0]
                g_species = S.get_leaves_by_name(name=g_species)[0]
                pureM[g] = g_species
                combinedM[g] = g_species
                smap[g] = g_species
                g.I = 1
                g.P = 1
            else:
                g.P = None  # init
                # computing pureM
                same_locus = [c for c in g.children if c.locus_id == g.locus_id]
                same_pureM = [pureM[c] for c in same_locus if pureM[c] is not None]
                if not same_locus or not same_pureM:
                    pureM[g] = None
                    warn("Detected an extinct lineage (node id %i); This may indicate corrupted data." % g.nid)
                    print g.get_ascii(attributes=['name', 'nid'])
                elif len(same_pureM) == 1:
                    pureM[g] = same_pureM[0]
                else:
                    pureM[g] = S.get_common_ancestor(same_pureM)
                # computing combinedM and I mapping
                all_pureM = [pureM[c] for c in g.children if pureM[c] is not None]
                if not all_pureM:
                    combinedM[g] = None
                    g.I = None
                elif len(all_pureM) == 1:
                    combinedM[g] = all_pureM[0]
                    g.I = combinedM[g].rank
                else:
                    combinedM[g] = S.get_common_ancestor(all_pureM)
                    g.I = combinedM[g].rank

        # P mapping
        for s in S.traverse():
            s.lastvisited = None

        leaves = [l for l in L]
        for i in range(1, d + 1):
            for lid, l1 in enumerate(leaves):
                if smap[l1].rank == i:
                    for l2 in leaves[(lid+1):]:
                        if smap[l2] == smap[l1]:
                            p = l1.get_common_ancestor(l2)
                            locus_set = [p.locus_id] + [c.locus_id for c in p.children]
                            if p.P is None and  {l1.locus_id, l2.locus_id}.issubset(locus_set):
                                p.P = i
                    smap[l1] = smap[l1].up

    def _label_source_nodes(self):
        """
        Labels the source nodes, i.e. the ends of cut edges, by adding a boolean 'source' attribute
        to each node.
        :return:
        """
        nodetypes = self._get_node_types()
        L = self.L
        for g in L.traverse():
            g.source = False
        L.source = True
        for i, g in enumerate(L.traverse(strategy="postorder")):
            if nodetypes[i] == "outside":
                if nodetypes[g.children[1].nid] == "root":
                    g.children[1].source = True
                else:
                    g.children[0].source = True

    def locus_tree(self):
        """
        Returns the gene tree with source nodes labelled.
        Source nodes represent the locations of evolutionary locus gain events in the gene tree.
        Formally, they are defined as the lower nodes of a cut edge.
        Each node in the returned tree contains a boolean `source` attribute.
        The nodes are also labelled with adjusted I and P mappings (computed based on "neighbouring" locus subtrees).
        :return: ete2.Tree object
        """
        if self.L:
            return self.L
        else:
            self.L = self.G.copy()
            self._label_source_nodes()
            self._assign_locus_ids()
            self._compute_adjusted_mappings()
            return self.L

    def forest(self):
        """
        Returns the decomposition as a forest of locus subtrees.
        :return: list
            A list of ete3.Tree objects.
        """
        if self.F:
            return self.F
        else:
            # G = self.G.copy()
            # for r in self.roots:
            #     self.F.append(G.search_nodes(nid=r)[0].detach())
            #     while len(G.children) == 1:  # moving down the root
            #         G = G.children[0]
            #     G.prune(G)
            # return self.F
            G = self.G.copy()
            self.F = [g.detach() for g in G.traverse() if g.nid in self.roots]
            self.F = [f.children[0] if len(f.children) == 1 else f for f in self.F]  # removing upper single nodes
            map(lambda x: x.prune(x), self.F)  # removing internal single nodes
            return self.F

    def roots(self):
        """
        Returns the list of subtree root IDs in postorder numbering.
        :return:
        """
        return self.roots

    def gene_tree(self):
        """
        Returns the gene tree, labelled with the raw I/P mapping values (i.e. without considering the locus
        structure in the tree).
        :return: ete2.Tree
        """
        if not self.G_labelled:
            compute_mappings(self.G, self.S)
        return self.G

    def show(self, layout=None, tree_style=None, palette=None):
        """
        Starts an interactive session to visualize the decomposition.
        :return: None
        """
        if not self.L:
            self.locus_tree()  # computes & stores the tree
        raise NotImplemented

    def render(self, fname, layout=None, tree_style=None, palette=None):
        """
        Renders the locus tree and writes the image to file.
        :param fname: str
            Output file path
        :param layout:
        :param tree_style:
        :param palette:
        :return: None
        """
        raise NotImplemented

    def write_forest(self, path=None):
        """
        Writes the decomposition forest in Newick format.
        The forest is returned as a string, optionally written to file.
        :param path: str
            Path to save the forest.
        :return: str
        """
        raise NotImplemented


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
            n.rank = 1
        else:
            n.rank = max(c.rank for c in n.children) + 1


def compute_mappings(G, S):
    """
    Computes both I and P mappings using Pawel's algorithm. Modifies G in situ.
    Leaf names in S are assumed to be unique.
    Leaf names in G are assumed to correspond to leaf names in S (no identifiers!)
    :param G: ete2.Tree
    :param S: ete2.Tree
    :return: None
    """
    try:
        d = S.rank
    except AttributeError:
        print("Species tree not initialized; Assign ranks before computing mappings.")
        raise
    else:
        S = S.copy()

    # I mapping
    for g in G.traverse(strategy="postorder"):
        if g.is_leaf():
            g_species = g.name
            g_species = S.get_leaves_by_name(name=g_species)[0]
            g.M = g_species
            g.I = 1
            g.P = 1
            g.smap = g_species
        else:
            g_species = S.get_common_ancestor(g.children[0].M, g.children[1].M)
            g.M = g_species
            g.I = g_species.rank
            g.P = None

    # P mapping
    for s in S.traverse():
        s.lastvisited = None

    for i in range(1, d+1):
        for g in G:
            if g.smap.rank == i:
                if g.smap.lastvisited is not None:
                    p = g.get_common_ancestor(g.smap.lastvisited)
                    if p.P is None:
                        p.P = i
                g.smap.lastvisited = g
                g.smap = g.smap.up


def assign_lineages(S):
    """
    Assigns the evolutionary lineages (as lists of nodes) to leaves of S.
    Works in situ, i.e. modifies S.
    The species tree needs to be labelled with taxonomic ranks by the assign_ranks() function.
    :param S: ete2.Tree
    :return: None
    """
    _, d = S.get_farthest_leaf()
    d = int(d) + 1
    for s in S:
        p = s
        lineage = [S]*d
        while not p.is_root():
            pp = p.up
            for i in range(p.rank-1, pp.rank-1):
                lineage[i] = p
            p = pp
        s.lineage = lineage


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
            if g.I == 1 or g.I > g.P:
                is_minimal[g] = True
    return [g for g in is_minimal if is_minimal[g]]


def cut_tree(g, S):
    """
    Finds an optimal cut under a node g with respect to the reference species tree S.
    An optimal cut is chosen from egdes adjacent to the children of g.
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
        Checks if rank r divides evolutionary lineages of two leaf sets
        (i.e. if their lineages have disjoint taxa in rank r-1).
        :param r: int
            The rank, an integer from 1 to height of S.
        :param leafset: set
            A set of leaf nodes (ete3.Tree objects).
        :return: bool
        """
        if r == 1:
            # This would mean that we check an internal node with I == 1,
            # so there is no embeddability.
            return False
        ancset1 = set([l.M.lineage[r-2] for l in leafset1])
        ancset2 = set([l.M.lineage[r-2] for l in leafset2])
        return ancset1.isdisjoint(ancset2)

    if g.I == g.P and g.I > 1:
        raise ValueError("The node is embeddable - nothing to cut!")
    elif g.is_leaf():
        raise ValueError("The node is a leaf - nothing to cut!")
    # Assigning candidate nodes for the root of the detached subtree
    dplc = g.children  # duplication-like candidates
    hgtc1 = dplc[0].children  # hgt-like candidates; either 0 or 2
    hgtc2 = dplc[1].children
    lc1 = len(hgtc1)
    lc2 = len(hgtc2)  # number of hgt-like candidates
    embeddable = [True, True] + [False]*(lc1 + lc2)  # g's children always satisfy embeddability condition
    loss_cost = [len(dplc[0].M) + len(dplc[1].M)]*2 + [0]*(lc1 + lc2)
    new_tree_size = [len(dplc[0]), len(dplc[1])] + [0]*(lc1 + lc2)

    # checking embeddability and loss cost
    leafset2 = set(dplc[1])
    for i in range(lc1):
        # if hgtc1 not empty, then it has two elements; 1-i is the other one
        leafset1 = set(hgtc1[1-i])
        newM = S.get_common_ancestor(dplc[1].M, hgtc1[1-i].M)  # updated M mapping after removal of hgtc1[i]
        newI = newM.rank
        embeddable[2+i] = check_division(leafset1, leafset2, newI)
        loss_cost[2+i] = len(hgtc1[i].M) + len(newM)  # turned out that investigating size of G is unneccessary
        new_tree_size[2+i] = len(hgtc1[i])

    leafset2 = set(dplc[0])
    for i in range(lc2):
        leafset1 = set(hgtc2[1-i])
        newM = S.get_common_ancestor(dplc[0].M, hgtc2[1-i].M)
        newI = newM.rank
        embeddable[2 + lc1 + i] = check_division(leafset1, leafset2, newI)
        loss_cost[2 + lc1 + i] = len(hgtc2[i].M) + len(newM)
        new_tree_size[2 + lc1 + i] = len(hgtc2[i])

    min_cost = min(loss_cost)
    good_cut = [loss_cost[i] == min_cost and embeddable[i] for i in range(2+lc1+lc2)]

    cut_id = 0
    for i in range(2 + lc1 + lc2):
        if good_cut[i] and new_tree_size[i] < new_tree_size[cut_id]:
            # We cut out the smallest tree.
            # It makes the forest slightly higher but lowers the number of non-admissible events.
            cut_id = i

    if cut_id < 2:
        newroot = dplc[cut_id]
    elif cut_id < 2 + lc1:
        newroot = hgtc1[cut_id - 2]
    elif cut_id < 2 + lc1 + lc2:
        newroot = hgtc2[cut_id - lc1 - 2]
    else:
        raise RuntimeError("Wrong rootnode id!")
    return newroot


def decompose(gene_tree, species_tree):
    """
    Performs the decomposition of the gene tree with respect to the species tree.
    Returns the forest of subtrees as a list.
    :param gene_tree: ete2.Tree object
    :param species_tree:  ete2.Tree object
        The reference species tree. The leaf names need to be unique.
    :return: list
        A list of ete2.Tree objects, each one corresponding to one locus subtree from
        the decomposition forest.
    """
    G = gene_tree.copy()
    S = species_tree.copy()

    # validate the data
    snames = set()
    for s in S:
        if s.name in snames:
            raise ValueError("Node names of the species tree are not unique!")
        else:
            snames.add(s.name)

    # names of leaves are stripped from the identifiers
    # original names of leaves in G are stored, to be assigned back at the end of the function
    rawnames = dict((g, g.name) for g in G)
    for g in G:
        g.name = g.name.split('_')[0]
        if g.name not in snames:
            raise ValueError("The gene %s does not correspond to any species!" % g.name)

    for i, g in enumerate(G.traverse(strategy="postorder")):
        g.id = i

    # initialize data
    assign_ranks(S)
    assign_lineages(S)

    # decompose
    improper = [G]  # list of minimal "improper" nodes
    roots = []  # decomposition forest
    while improper:
        compute_mappings(G, S)
        improper = minimal_nodes(G)
        subtrees = [cut_tree(g, S).detach() for g in improper]
        roots += [s.id for s in subtrees]
        while len(G.children) == 1:
            G = G.children[0]  # moving down the root, because pruning preserves it
        G.prune(G)  # removing nodes with single children; whether the node or it's child is removed doesn't matter
                    # because their parents are not improper anymore
    roots.append(G.id)
    return Decomposition(gene_tree.copy(), S, roots)


if __name__=="__main__":
    from time import time
    S = ete2.Tree('(((a, b), (c, d)), (e, (f, g)));')
    G = ete2.Tree('((((a_1, c_1), b_1), (e_1, f_1)), (c_2, ((g_1, f_2), (d_1, a_2))));')
    #S = ete2.Tree('(a, b);')
    #G = ete2.Tree('((a, b), (a, b));')
    #S = ete2.Tree('(((a, b), (c, d)), d);')
    #G = ete2.Tree('((((a, b), (c, d)), (a, b)), d);')
    #S = ete2.Tree('/home/ciach/Projects/Tree_node_classification/ScoringSystem/TransferOnly/Loss00/S00/species_tree.labelled.tree', format=8)
    #G = ete2.Tree('/home/ciach/Projects/Tree_node_classification/ScoringSystem/TransferOnly/Loss00/S00/gene_tree2.simdec.tree')
    #G = ete2.Tree('/home/ciach/Projects/Tree_node_classification/ScoringSystem/BothEvents/Loss09/S07/gene_tree8.simdec.tree', format=9)
    #S = ete2.Tree('/home/ciach/Projects/Tree_node_classification/ScoringSystem/BothEvents/Loss09/S07/species_tree.labelled.tree', format=8)
    # S = ete2.Tree('((v, s), (x, w));')
    # G = ete2.Tree('(((s, v), x), (w, w));')
    #print S
    #print G

    #S = ete2.Tree('/home/ciach/Projects/TreeDecomposition/Simulations/DTL_stochastic/S04/species_tree.labelled.tree', format=8)
    #G = ete2.Tree('/home/ciach/Projects/TreeDecomposition/Simulations/DTL_stochastic/S04/gene_tree6.labelled.tree', format=9)

    #G = G.get_common_ancestor(G&"S323_1", G&"S313_1")
    #S = S.get_common_ancestor([S&g.name.split('_')[0] for g in G])

    s = time()
    D = decompose(G, S)
    e = time()
    print "Decomposed in %.02f seconds" % (e - s)
    s = time()
    F = D.forest()
    e = time()
    print "Forest obtained in %.02f seconds" % (e - s)
    s = time()
    L = D.locus_tree()
    e = time()
    print "Locus tree obtained in %.02f seconds" % (e - s)
