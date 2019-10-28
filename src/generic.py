"""
Methods to handle decomposition results regardless of algorithm used
"""
from warnings import warn
import re

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
            
class Decomposition(object):
    def __init__(self, gene_tree, species_tree, roots):
        """
        A locus decomposition of a gene tree with respect to a given species tree.
        :param gene_tree: ete2.Tree
            A gene tree. Each leaf name needs to correspond exactly to one leaf name of the species tree.
        :param species_tree: ete2.Tree
            A ranked species tree. Each node nees to have a 'rank' attribute.
            If the supplied species tree is not ranked, it is possible to assign artificial
            ranks after initializing the object by running
            the method assign_topological_ranks() from this class.
            The ranks need to be assigned prior to most operations on the decomposition.
        :param roots: list
            A list of identifiers (in post-order numbering) of the roots of trees from the decomposition forest.
        """
        self.G = gene_tree.copy()
        for i, g in enumerate(self.G.traverse(strategy="postorder")):
            g.nid = i  # used to identify e.g. forest roots
        self.G_labelled = False  # whether G has been labelled with raw I/P mappings
        self.S = species_tree.copy()
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
            matching_species = []
            for s in self.S:
                if re.match(s.name, l.name):
                    matching_species.append(s)
            matching_nb = len(matching_species)
            if matching_nb == 0:
                raise ValueError('Gene tree leaf %s does not correspond to any species!' % original_gene_names[i])
            elif matching_nb > 1:
                raise ValueError('Ambiguous species mapping for gene tree leaf %s' % (original_gene_names[i], matching_names[0], matching_names[1]))
            else:
                l.species = matching_species[0]
##        for l in self.G:
##            l.species = self.S.get_leaves_by_name(l.name.split('_')[0])
##            assert len(l.species) == 1, "Species names are not unique!"
##            l.species = l.species[0]


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
                    print(g.get_ascii(attributes=['name', 'nid']))
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

        # def _score(nid1, nid2):
        #     """
        #     Score a path joining two roots from the forest.
        #     The roots are identified by their postorder IDs: nid1, nid2.
        #     :return: float
        #     """
        #     return 1 if {nid1, nid2} in pairs else 0
        def _score(node, nid1, nid2):
            """
            Score a path joining two roots from the forest.
            The roots are identified by their postorder IDs: nid1, nid2.
            :return: float
            """
            r1 = L.search_nodes(nid = nid1)[0]
            r2 = L.search_nodes(nid = nid2)[0]
            sp1 = [l.species for l in r1]
            sp2 = [l.species for l in r2]
            I1 = sp1[0].get_common_ancestor(sp1).rank
            I2 = sp2[0].get_common_ancestor(sp2).rank
            I = sp1[0].get_common_ancestor(sp1+sp2).rank
            phi = I - (I1 + I2)/2.
            phi /= self.S.rank
            # _, depth = node.get_farthest_leaf(topology_only=True)
            # depth += 1
            return (1 if {nid1, nid2} in pairs else 0, 1-phi)

        def _psum(*tuples):
            return tuple(map(sum, zip(*tuples)))

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
                g.pcosts = [(0., 0.)]  # partial costs
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
                    m = max(_psum(c2.pcosts[i2],  _score(g, r1, r2)) for i2, r2 in enumerate(c2.roots))
                    g.pcosts.append(_psum(c1.pcosts[i1], m))
                    g.optimals.append(c1.optimals[i1] * sum(c2.optimals[i2] for
                                                            i2, r2 in enumerate(c2.roots) if
                                                            _psum(c2.pcosts[i2], _score(g, r1, r2)) == m))
                for i2, r2 in enumerate(c2.roots):
                    m = max(_psum(c1.pcosts[i1], _score(g, r2, r1)) for i1, r1 in enumerate(c1.roots))
                    g.pcosts.append(_psum(c2.pcosts[i2], m))
                    g.optimals.append(c2.optimals[i2] * sum(c1.optimals[i1] for
                                                            i1, r1 in enumerate(c1.roots) if
                                                            _psum(c1.pcosts[i1], _score(g, r2, r1)) == m))
        self.optimal_score = (0, 0)
        # backtracking
        for g in L.traverse(strategy="preorder"):
            if g.inspect and (g.is_root() or not g.up.inspect):
                # g is a root of a junction node cluster
                g.source = True
                if g.is_leaf() and not g.roots == [g.nid]:
                    raise ValueError("Leaf not a subtree root")
                if g.nid in self.roots:
                    # g is a root of a locus subtree
                    continue
                else:
                    c1, c2 = g.children
                    full_cost = (-1, -1)
                    root1 = -1
                    root2 = -1
                    for i1, r1 in enumerate(c1.roots):
                        for i2, r2 in enumerate(c2.roots):
                            ncost = _psum(c1.pcosts[i1], c2.pcosts[i2], _score(g, r1, r2))
                            if ncost > full_cost:
                                full_cost = ncost
                                root1 = r1
                                root2 = r2
                    assert root1 >= 0 and root2 >= 0
                    self.optimal_score = _psum(full_cost, self.optimal_score)
                    g.cluster_opts = sum(c1.optimals[i1] * c2.optimals[i2] for
                                         i1, r1 in enumerate(c1.roots) for i2, r2 in enumerate(c2.roots) if
                                         _psum(c1.pcosts[i1], c2.pcosts[i2], _score(g, r1, r2)) == full_cost)
                    c2.source = True  # arbitrary choice of source below a cluster root
                    g.root = root1  # joining g with a root by a path with no cuts
                    c1.root = root1
                    c2.root = root2
            elif g.inspect:
                # "internal" node of a junction cluster
                if g.is_leaf() and not g.roots == [g.nid]:
                    raise ValueError("Leaf not a subtree root")
                if g.nid in self.roots:
                    continue
                c1, c2 = g.children
                if g.root in c1.roots:
                    c2.source = True
                    c1.root = g.root
                    source_child = c2
                elif g.root in c2.roots:
                    c1.source = True
                    c2.root = g.root
                    source_child = c1
                else:
                    raise ValueError("Node joined to non-reachable subtree")
                max_cost = (-1, -1)
                new_root = -1
                for i, r in enumerate(source_child.roots):
                    ncost = _psum(source_child.pcosts[i], _score(g, g.root, r))
                    if ncost > max_cost:
                        new_root = r
                        max_cost = ncost
                assert new_root >= 0
                source_child.root = new_root
        self.number_of_optimal_solutions = 1
        for g in L.traverse():
            self.number_of_optimal_solutions *= g.cluster_opts 


    def _assert_sources(self):
        internal_nids = [n.nid for f in self.F for n in f.iter_descendants()]
        root_nids = [f.nid for f in self.F]
        outside_nids = [n.nid for n in self.G.traverse() if n.nid not in internal_nids + root_nids]
        source_nids = [n.nid for n in self.L.traverse() if n.source]
        junction_nids = [n.nid for n in self.L.traverse() if any(c.source for c in n.children)]
        assert junction_nids == outside_nids

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
