import ete2
from warnings import warn


class Decomposition(object):
    def __init__(self, gene_tree, species_tree, roots):
        """
        A locus decomposition of a gene tree with respect to a given species tree.
        :param gene_tree: ete2.Tree
        :param species_tree: ete2.Tree
            A ranked species tree. Each node nees to have a 'rank' attribute.
            If the supplied species tree is not ranked, it is possible to assign artificial
            ranks after initializing the object by running
            the method assign_topological_ranks() from this class.
            The ranks need to be assigned prior to most operations on the decomposition.
        :param roots: list
            A list of identifiers (in post-order numbering) of the roots of trees from the decomposition forest.
        """
        self.G = gene_tree
        for i, g in enumerate(self.G.traverse(strategy="postorder")):
            g.nid = i  # used to identify e.g. forest roots
        self.G_labelled = False  # whether G has been labelled with raw I/P mappings
        self.S = species_tree
        self.roots = sorted(roots)
        self.F = []  # forest; self.F[i] corresponds to self.roots[i]
        self.L = None  # locus tree
        self.colorized = False  # whether the locus tree has been colorized for rendering
        self._map_gene_to_species()
        # self._assign_locus_ids()

    def assign_topological_ranks(self):
        """
        Assign ranks to the nodes of the species tree.
        A topological rank is equal to the depth of the node.
        :return: None
        """
        for n in self.S.traverse(strategy="postorder"):
            if n.is_leaf():
                n.rank = 0
            else:
                n.rank = max(c.rank for c in n.children) + 1

    def _map_gene_to_species(self):
        for l in self.G:
            l.species = self.S.get_leaves_by_name(l.name.split('_')[0])
            assert len(l.species) == 1;
            "Species names are not unique!"
            l.species = l.species[0]

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
                g.I = 0
                g.P = 0
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
        for i in range(0, d + 1):
            for lid, l1 in enumerate(leaves):
                if smap[l1].rank == i:
                    for l2 in leaves[(lid + 1):]:
                        if smap[l2] == smap[l1]:
                            p = l1.get_common_ancestor(l2)
                            locus_set = [p.locus_id] + [c.locus_id for c in p.children]
                            if p.P is None and {l1.locus_id, l2.locus_id}.issubset(locus_set):
                                p.P = i
                    smap[l1] = smap[l1].up

    def _lift_roots(self):
        """
        Labels the source nodes, i.e. the ends of cut edges, by adding a boolean 'source' attribute
        to each node of the locus tree.
        Maximizes the number of junction nodes with adjusted P mapping value = 0.
        :return:
        """

        def _score(nid1, nid2):
            """
            Score a path joining two roots from the forest.
            The roots are identified by their postorder IDs: nid1, nid2.
            :return: float
            """
            return 1 if {nid1, nid2} in pairs else 0

        L = self.L
        if not L:
            raise ValueError("Locus tree not initialized. Use the locus_tree() method.")
        # identifying forests sharing a species
        forest = self.forest()
        labels = [{l.species for l in f} for f in forest]
        pairs = [{forest[i].nid, forest[j].nid} for i in range(len(forest))
                 for j in range(i + 1, len(forest)) if not labels[i].isdisjoint(labels[j])]

        # computing partial costs
        for g in L.traverse():
            g.cluster_opts = 1  # nb of optimal solutions in cluster
            g.source = False
            if g.nid in self.roots:
                g.inspect = True  # whether to compute cost for g
                g.roots = [g.nid]  # postorder IDs of forest roots visible from g
                g.pcosts = [0]  # partial costs
                g.root = g.nid  # index of root joined with g by a path with no cuts
                g.optimals = [1]  # number of optimal solutions
            else:
                g.roots = []
                g.pcosts = []
                g.inspect = False
                g.root = -1
                g.optimals = []
        for g in L.traverse(strategy="postorder"):
            if g.is_leaf():
                continue
            c1, c2 = g.children
            if c1.inspect and c2.inspect:
                g.inspect = True
                g.roots = c1.roots + c2.roots
                for i1, r1 in enumerate(c1.roots):
                    m = max(c2.pcosts[i2] + _score(r1, r2) for i2, r2 in enumerate(c2.roots))
                    g.pcosts.append(c1.pcosts[i1] + m)
                    g.optimals.append(c1.optimals[i1] * sum(c2.optimals[i2] for
                                                            i2, r2 in enumerate(c2.roots) if
                                                            c2.pcosts[i2] + _score(r1, r2) == m))
                for i2, r2 in enumerate(c2.roots):
                    m = max(c1.pcosts[i1] + _score(r1, r2) for i1, r1 in enumerate(c1.roots))
                    g.pcosts.append(c2.pcosts[i2] + m)
                    g.optimals.append(c2.optimals[i2] * sum(c1.optimals[i1] for
                                                            i1, r1 in enumerate(c1.roots) if
                                                            c1.pcosts[i1] + _score(r1, r2) == m))
        # computing full costs & backtracking
        for g in L.traverse():
            if g.inspect and (g.is_root() or not g.up.inspect):
                # g is a root of a junction node cluster
                g.source = True
                if g.is_leaf():
                    continue
                elif g.roots == [g.nid]:
                    continue  # g is a root of a subtree
                c1, c2 = g.children
                full_cost = -1
                root1 = -1
                root2 = -1
                for i1, r1 in enumerate(c1.roots):
                    for i2, r2 in enumerate(c2.roots):
                        ncost = c1.pcosts[i1] + c2.pcosts[i2] + _score(r1, r2)
                        if ncost > full_cost:
                            full_cost = ncost
                            root1 = r1
                            root2 = r2
                g.cluster_opts = sum(c1.optimals[i1] * c2.optimals[i2] for
                                     i1, r1 in enumerate(c1.roots) for i2, r2 in enumerate(c2.roots) if
                                     c1.pcosts[i1] + c2.pcosts[i2] + _score(r1, r2) == full_cost)
                g.cluster_opts *= 2  # to account for no constraint for the choice of cut below the root
                r1 = L.search_nodes(nid=root1)[0]  # note that order of trees in forest corresponds to self.roots
                r2 = L.search_nodes(nid=root2)[0]  # this is time consuming, probably should be implemented faster
                assert r1.root == root1
                assert r2.root == root2
                while r1 != g:
                    s = r1.get_sisters()[0]
                    s.source = True
                    r1.root = root1
                    r1 = r1.up
                while r2 != g:
                    s = r2.get_sisters()[0]
                    s.source = True
                    r2.root = root2
                    r2 = r2.up
                c1.source = False  # both children labelled as source, need to erase one labelling
                g.root = root1  # joining g with a root by a path with no cuts
                # traversing the junction cluster below g:
                for v in g.traverse(strategy="preorder",
                                    is_leaf_fn=lambda x: x.nid == x.root):
                    if v.nid == v.root: continue  # root of a forest tree
                    assert not all(c.source for c in v.children), "Both children are source nodes!"
                    c = v.children[0] if v.children[0].source else v.children[1]
                    mcost = -1
                    croot = -1
                    for i, r in enumerate(c.roots):
                        ncost = c.pcosts[i] + _score(v.root, r)
                        if ncost > mcost:
                            mcost = ncost
                            croot = r
                    r = L.search_nodes(nid=croot)[0]
                    assert r.root == croot
                    c.root = croot
                    while (r != c):
                        s = r.get_sisters()[0]
                        s.source = True
                        r.root = croot
                        r = r.up
        self.number_of_optimal_solutions = reduce(lambda x, y: x*y, [g.cluster_opts for g in L.traverse()])

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
            self._lift_roots()
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
            G = self.G.copy()
            for r in self.roots:
                self.F.append(G.search_nodes(nid=r)[0].detach())
                while len(G.children) == 1:  # pushing down the root
                    G = G.children[0]
                G.prune(G)
            return self.F

    def gene_tree(self):
        """
        Returns the gene tree, labelled with the raw I/P mapping values (i.e. without considering the locus
        structure in the tree).
        :return: ete2.Tree
        """
        if not self.G_labelled:
            compute_mappings(self.G, self.S)
        return self.G

    def _colorize(self, palette):
        """
        Assigns faces and colours to the locus trees for pretty rendering.
        :param palette: list
            List of strings representing colours in hexadecimal format.
        :return:
        """
        from ete2 import NodeStyle, faces
        if not self.L:
            self.locus_tree()  # computes & stores the tree
        ncol = len(palette)
        iFace = faces.AttrFace("I", fsize=8, text_suffix='/')
        pFace = faces.AttrFace("P", fsize=8)
        # idFace = faces.AttrFace("id", fsize=8)
        # suppFace = faces.AttrFace("support", text_suffix=" ", formatter="%.2f", fsize=8)
        coloured = dict((i, False) for i, g in enumerate(self.L.traverse(strategy="postorder")))
        current_colour = -1
        for g in self.L.traverse(strategy="postorder"):
            if not g.is_leaf():
                # g.add_face(suppFace, position="branch-bottom", column=-2)
                g.add_face(iFace, position="branch-top", column=-1)
                g.add_face(pFace, position="branch-top", column=0)
            if g.source:
                current_colour += 1
                current_colour %= ncol
                style = NodeStyle()
                style['vt_line_color'] = palette[current_colour]
                style['hz_line_color'] = palette[current_colour]
                style['size'] = 0
                style['fgcolor'] = '#000000'
                style["vt_line_width"] = 2
                style["hz_line_width"] = 2
                for gg in g.traverse():
                    if not coloured[gg.nid]:
                        gg.set_style(style)
                        coloured[gg.nid] = True

    def show(self, tree_style=None, palette=None):
        """
        Starts an interactive session to visualize the decomposition.
        :return: None
        """
        if not palette:
            palette = ['#1F77B4', '#AEC7E8', '#FF7F0E',
                       '#FFBB78', '#2CA02C', '#98DF8A',
                       '#D62728', '#FF9896', '#9467BD',
                       '#C5B0D5', '#8C564B', '#C49C94',
                       '#E377C2', '#F7B6D2', '#7F7F7F',
                       '#C7C7C7', '#BCBD22', '#DBDB8D',
                       '#17BECF', '#9EDAE5']
        if not self.colorized:
            self._colorize(palette)
            self.colorized = True
        if not tree_style:
            from ete2 import TreeStyle
            tree_style = TreeStyle()
            # tstyle.show_leaf_name = False
            tree_style.scale = 28
            tree_style.branch_vertical_margin = 6
            tree_style.show_branch_length = False

            # tstyle.show_branch_support = True
            tree_style.show_scale = False
        self.L.convert_to_ultrametric()
        self.L.show(tree_style=tree_style)

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
        if not palette:
            palette = ['#1F77B4', '#AEC7E8', '#FF7F0E',
                       '#FFBB78', '#2CA02C', '#98DF8A',
                       '#D62728', '#FF9896', '#9467BD',
                       '#C5B0D5', '#8C564B', '#C49C94',
                       '#E377C2', '#F7B6D2', '#7F7F7F',
                       '#C7C7C7', '#BCBD22', '#DBDB8D',
                       '#17BECF', '#9EDAE5']
        if not self.colorized:
            self._colorize(palette)
            self.colorized = True
        if not tree_style:
            from ete2 import TreeStyle
            tree_style = TreeStyle()  # imported during colorizing tree
            # tstyle.show_leaf_name = False
            tree_style.scale = 28
            tree_style.branch_vertical_margin = 6
            tree_style.show_branch_length = False

            # tstyle.show_branch_support = True
            tree_style.show_scale = False
        self.L.convert_to_ultrametric()
        self.L.render(file_name=fname, tree_style=tree_style)

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
            n.rank = 0
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
        lineage = [S] * d
        while not p.is_root():
            pp = p.up
            for i in range(p.rank - 1, pp.rank - 1):
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
            if g.I == 0 or g.I > g.P:
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
            The rank, an integer from 0 to height of S.
        :param leafset: set
            A set of leaf nodes (ete3.Tree objects).
        :return: bool
        """
        if r == 0:
            # This would mean that we check an internal node with I == 0 after pruning,
            # so there is no embeddability.
            return False
        ancset1 = set([l.M.lineage[r - 1] for l in leafset1])
        ancset2 = set([l.M.lineage[r - 1] for l in leafset2])
        return ancset1.isdisjoint(ancset2)

    if g.I == g.P and g.I > 0:
        raise ValueError("The node is embeddable - nothing to cut!")
    elif g.is_leaf():
        raise ValueError("The node is a leaf - nothing to cut!")
    # Assigning candidate nodes for the root of the detached subtree
    dplc = g.children  # duplication-like candidates
    hgtc1 = dplc[0].children  # hgt-like candidates; either 0 or 2
    hgtc2 = dplc[1].children
    lc1 = len(hgtc1)
    lc2 = len(hgtc2)  # number of hgt-like candidates
    embeddable = [True, True] + [False] * (lc1 + lc2)  # g's children always satisfy embeddability condition
    loss_cost = [len(dplc[0].M) + len(dplc[1].M)] * 2 + [0] * (lc1 + lc2)
    new_tree_size = [len(dplc[0]), len(dplc[1])] + [0] * (lc1 + lc2)

    # checking embeddability and loss cost
    leafset2 = set(dplc[1])
    for i in range(lc1):
        # if hgtc1 not empty, then it has two elements; 1-i is the other one
        leafset1 = set(hgtc1[1 - i])
        newM = dplc[1].M.get_common_ancestor(hgtc1[1 - i].M)  # updated M mapping after removal of hgtc1[i]
        newI = newM.rank
        embeddable[2 + i] = check_division(leafset1, leafset2, newI)
        loss_cost[2 + i] = len(hgtc1[i].M) + len(newM)  # turned out that investigating size of G is unneccessary
        new_tree_size[2 + i] = len(hgtc1[i])

    leafset2 = set(dplc[0])
    for i in range(lc2):
        leafset1 = set(hgtc2[1 - i])
        newM = dplc[0].M.get_common_ancestor(hgtc2[1 - i].M)
        newI = newM.rank
        embeddable[2 + lc1 + i] = check_division(leafset1, leafset2, newI)
        loss_cost[2 + lc1 + i] = len(hgtc2[i].M) + len(newM)
        new_tree_size[2 + lc1 + i] = len(hgtc2[i])

    # not very efficient, but oh well
    min_cost = min([loss_cost[i] for i in range(2 + lc1 + lc2) if embeddable[i]])
    # embeddable & min cost:
    good_cut = [loss_cost[i] == min_cost and embeddable[i] for i in range(2 + lc1 + lc2)]
    target_size = min([new_tree_size[i] for i in range(2 + lc1 + lc2) if good_cut[i]])
    # embeddable, min cost & min tree:
    cut_id = [i for i in range(2 + lc1 + lc2) if good_cut[i] and new_tree_size[i] == target_size][0]

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
        g.nid = i

    # initialize data
    assign_ranks(S)
    assign_lineages(S)

    # decompose
    improper = [G]  # list of minimal non-embeddable nodes
    roots = []  # decomposition forest
    while improper:
        compute_mappings(G, S)
        improper = minimal_nodes(G)
        subtrees = [cut_tree(g, S) for g in improper]
        roots += [s.nid for s in subtrees]
        map(lambda x: x.detach(), subtrees)
        # Pruning G. This is time costly, and ideally should be get rid of,
        # but is required in the current implementation.
        while len(G.children) == 1:
            G = G.children[0]  # moving down the root, because pruning preserves it
        G.prune(G)  # removing nodes with single children
    roots.append(G.nid)
    return Decomposition(gene_tree, S, roots)


if __name__ == "__main__":
    from time import time

    # S = ete2.Tree('(((a, b), (c, d)), (e, (f, g)));')
    # G = ete2.Tree('((((a_1, c_1), b_1), (e_1, f_1)), (c_2, ((g_1, f_2), (d_1, a_2))));')
    # S = ete2.Tree('(a, b);')
    # G = ete2.Tree('((a, b), (a, b));')
    # S = ete2.Tree('(((a, b), (c, d)), d);')
    # G = ete2.Tree('((((a, b), (c, d)), (a, b)), d);')
    # S = ete2.Tree("(a, b);")
    # G = ete2.Tree("(((a, a), (a, a)), (b, b));")
    # S = ete2.Tree('/home/ciach/Projects/Tree_node_classification/ScoringSystem/TransferOnly/Loss00/S00/species_tree.labelled.tree', format=8)
    # G = ete2.Tree('/home/ciach/Projects/Tree_node_classification/ScoringSystem/TransferOnly/Loss00/S00/gene_tree2.simdec.tree')
    # G = ete2.Tree('/home/ciach/Projects/Tree_node_classification/ScoringSystem/BothEvents/Loss09/S07/gene_tree8.simdec.tree', format=9)
    # S = ete2.Tree('/home/ciach/Projects/Tree_node_classification/ScoringSystem/BothEvents/Loss09/S07/species_tree.labelled.tree', format=8)
    # S = ete2.Tree('((v, s), (x, w));')
    # G = ete2.Tree('(((s, v), x), (w, w));')
    # print S
    # print G

    # S = ete2.Tree('/home/ciach/Projects/TreeDecomposition/Simulations/DTL_stochastic/S04/species_tree.labelled.tree', format=8)
    # G = ete2.Tree('/home/ciach/Projects/TreeDecomposition/Simulations/DTL_stochastic/S04/gene_tree6.labelled.tree', format=9)

    # G = G.get_common_ancestor(G&"S228_1", G&"S290_1")
    # S = S.get_common_ancestor([S&g.name.split('_')[0] for g in G])

    # G = ete2.Tree('((((((((S21_1:1,(S22_3:1,(S22_9:1,S22_10:1)1:1)1:1)1:1,((S24_1:1,(S25_1:1,S26_1:1)1:1)1:1,((S29_9:1,(S35_14:1,S6_4:1)1:1)1:1,(S0_5:1,S0_6:1)1:1)1:1)1:1)1:1,S32_1:1)1:1,(((((S21_7:1,(S21_11:1,(S22_15:1,S38_7:1)1:1)1:1)1:1,S21_4:1)1:1,(S22_5:1,S22_6:1)1:1)1:1,(((S25_3:1,S26_3:1)1:1,S29_6:1)1:1,(S25_2:1,(S49_3:1,S26_4:1)1:1)1:1)1:1)1:1,(S32_3:1,S32_4:1)1:1)1:1)1:1,((S34_4:1,S35_2:1)1:1,S35_1:1)1:1)1:1,(((S38_3:1,((S63_3:1,(S64_3:1,S65_3:1)1:1)1:1,(S24_3:1,S43_4:1)1:1)1:1)1:1,((S39_2:1,S39_3:1)1:1,((S40_1:1,S41_1:1)1:1,(S43_3:1,S44_1:1)1:1)1:1)1:1)1:1,(S49_1:1,((S50_2:1,S50_3:1)1:1,(S51_2:1,S29_2:1)1:1)1:1)1:1)1:1)1:1,(((S56_1:1,S57_1:1)1:1,((S60_2:1,S60_3:1)1:1,S34_2:1)1:1)1:1,((S56_2:1,(S29_1:1,S57_3:1)1:1)1:1,((S63_1:1,(S64_1:1,S65_1:1)1:1)1:1,(S63_2:1,(S64_2:1,S65_2:1)1:1)1:1)1:1)1:1)1:1)1:1,(S19_2:1,(S19_3:1,(((((((S21_9:1,S21_10:1)1:1,(S35_10:1,(((S25_5:1,S26_9:1)1:1,S22_17:1)1:1,(S22_18:1,S22_19:1)1:1)1:1)1:1)1:1,(((S43_6:1,S24_5:1)1:1,(S25_4:1,(S26_7:1,S26_8:1)1:1)1:1)1:1,S29_7:1)1:1)1:1,(S32_6:1,S32_7:1)1:1)1:1,(((S21_6:1,(S2_1:1,S22_12:1)1:1)1:1,(S29_8:1,(S56_4:1,(S57_8:1,S57_9:1)1:1)1:1)1:1)1:1,((S35_11:1,S35_12:1)1:1,(S61_1:1,S35_13:1)1:1)1:1)1:1)1:1,((S49_2:1,(S60_4:1,(S50_4:1,S51_3:1)1:1)1:1)1:1,(S39_4:1,((S40_2:1,(S41_4:1,S26_6:1)1:1)1:1,((S21_12:1,S43_7:1)1:1,S44_2:1)1:1)1:1)1:1)1:1)1:1,(((((S35_7:1,S35_8:1)1:1,S0_3:1)1:1,((S1_2:1,(S2_3:1,(S2_8:1,S2_9:1)1:1)1:1)1:1,S4_1:1)1:1)1:1,(((((S1_3:1,(S2_10:1,S1_6:1)1:1)1:1,(S34_6:1,S61_2:1)1:1)1:1,(((S0_7:1,S1_7:1)1:1,(S2_11:1,S60_5:1)1:1)1:1,S4_3:1)1:1)1:1,(((S1_5:1,(S39_5:1,(S2_15:1,S2_16:1)1:1)1:1)1:1,(S4_6:1,S4_7:1)1:1)1:1,S6_3:1)1:1)1:1,((S41_2:1,S8_1:1)1:1,(((S11_5:1,S11_6:1)1:1,(S12_2:1,S13_2:1)1:1)1:1,(S11_4:1,(((S12_4:1,S13_4:1)1:1,(S12_5:1,S13_5:1)1:1)1:1,(S12_3:1,S13_3:1)1:1)1:1)1:1)1:1)1:1)1:1)1:1,((((S11_2:1,S1_1:1)1:1,(S12_1:1,S13_1:1)1:1)1:1,(S56_3:1,(S57_5:1,S57_6:1)1:1)1:1)1:1,S19_4:1)1:1)1:1)1:1)1:1)1:1);')
    # S = ete2.Tree('(((S0:1,((((S1:1,S2:1)1:1,S4:1)1:1,S6:1)1:1,((S8:1,S9:1)1:1,(S11:1,(S12:1,S13:1)1:1)1:1)1:1)1:1)1:1,S19:1)1:1,(((((((S21:1,S22:1)1:1,((S24:1,(S25:1,S26:1)1:1)1:1,S29:1)1:1)1:1,S32:1)1:1,(S34:1,S35:1)1:1)1:1,((S38:1,(S39:1,((S40:1,S41:1)1:1,(S43:1,S44:1)1:1)1:1)1:1)1:1,(S49:1,(S50:1,S51:1)1:1)1:1)1:1)1:1,(S56:1,S57:1)1:1)1:1,(((S60:1,S61:1)1:1,(S63:1,(S64:1,S65:1)1:1)1:1)1:1,S69:1)1:1)1:1);')
    # roots = [1, 2, 14, 17, 24, 25, 27, 34, 40, 44, 50, 57, 61, 66, 67, 72, 85, 89, 98, 105, 113, 124, 125, 128, 131, 134, 140, 144, 152, 157, 163, 166, 168, 169, 172, 179, 188, 191, 200, 206, 207, 213, 217, 220, 221, 223, 224, 228, 235, 236, 241, 244, 248, 251, 257, 261, 274, 276, 281, 283, 288]

    G = ete2.Tree('(((a, b), ((a, c), (c, d))), ((b, d), f));')
    S = ete2.Tree('(a, b, c, d, e, f);')
    roots = [0, 1, 5, 8, 13, 14]

    TD = Decomposition(G, S, roots)
    TD.assign_topological_ranks()
    TD.forest()
    TD.locus_tree()

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

    print L.get_ascii(attributes=['name', 'nid', 'source', 'type', 'I', 'P'])

    # for i in range(1000):
    #     print i
    #     G = ete2.Tree()
    #     G.populate(100, names_library=[s.name for s in S], reuse_names=True)
    #     D = decompose(G, S)
    #     F = D.forest()
    #     L = D.locus_tree()
    #     lnms = [set([l.name for l in L.search_nodes(locus_id=i) if l.is_leaf()]) for i in set(l.locus_id for l in L)]
    #     fnms = [set([l.name for l in f]) for f in F]
    #
    #     assert all([l in fnms for l in lnms])
    #     assert all([l in lnms for l in fnms])
