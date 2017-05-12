from operator import itemgetter
from itertools import groupby
import getopt, sys, ete2
from warnings import warn
import string
from matplotlib import pyplot as plt
import numpy as np
from time import time

# Data obtaining & curating
def check_data(lineages, max_unknown=3):
    """
    Checks if the lineages contain at most max_unknown UNKNOWN values.
    :param lineages:
    :return:
    """
    unknown_nbs = map(lambda x: sum(xx in {"UNKNOWN", '0'} for xx in x), lineages.itervalues())
    if max(unknown_nbs) > max_unknown:
        raise ValueError("""Classifying nodes of a lineage with more than %i UNKNOWN values is unreliable.
                            Please remove the genes with ambiguous lineages or change the missing data constraint.""" % max_unknown)


def curate_lineages(lineages):
    """
    Curates the evolutionary lineages: sets each UNKNOWN value to be an ancestor of the nearest lower rank.
    :param lineages: dict
        A dictionary of evolutionary lineages, indexed by TaxIDs.
    :return: dict
        A dictionary of curated lineages.
    """
    lineages = lineages.copy()
    for taxid in lineages:
        lge = lineages[taxid]
        for i, tx in enumerate(lge):
            if tx in {"UNKNOWN", '0'}:
                if i == 0:
                    lge[i] = tx + '_anc'
                else:
                    lge[i] = lge[i-1] + '_anc'
    lengths = map(len, lineages.itervalues())
    if len(set(lengths)) != 1:
        raise ValueError("The supplied lineages have unequal lengths.")
    return lineages


def annotate_tree(gene_tree, lineage_dict):
    """Annotates a gene tree, adding a lineage component to each leaf.
    Lineage component is a list containing a sequence of taxa - from lowest to highest - that
    represent the taxonomic classification of the leaf.

    Parameters
    ----------
    gene_tree: ete2.Tree object
        The gene tree to be annotated, with species' taxids as leaf names, possibly with a number after an underscore.
    lineages: list
        A dictionary of evolutionary lineages, indexed by TaxIDs.
    ranks: list
        List of evolutionary ranks (e.g. ['species', 'genus', 'kingdom']). The ranks should correspond to the
        keys of lineage dictionaries. Their purpose is to provide order of the keys.

    Value
    ----------
        ete2.Tree object equal to gene_tree with 'lineage' attribute added to each leaf.
        """
    gene_tree = gene_tree.copy()
    for leaf in gene_tree:
        taxid = leaf.name
        # stripping the number
        if '_' in taxid:
            taxid = taxid[:taxid.index('_')]
        leaf.lineage = lineage_dict[taxid]
    return gene_tree


def check_annotation(tree):
    """Checks the validity of gene tree taxonomic labelling.
    If a leaf has no 'lineage' attribute or there is no 'Biont' taxon in the lineage,
    then a RuntimeError is raised.
    :param tree: ete2.TreeNode object
        Annotated gene tree (i.e. each leaf has to contain a 'lineage' attribute that is
        a list of strings).
    :return: None
    """
    T = tree.copy()
    for leaf in T:
        if not hasattr(leaf, "lineage"):
            raise ValueError("Improperly labelled gene tree - no lineage for leaf %s." % leaf.name)

    # for leaf in T:
    #     # '131567' is defined as the biont taxID in shared.py and added manually to the lineage.
    #     # However, note that NCBI is sometimes unpredictable, so this taxID may be assigned to a
    #     # different taxon in the NCBI database in the future, so this is a possible source of bugs.
    #     if 'Biont' not in leaf.lineage and '131567' not in leaf.lineage:
    #         raise RuntimeError("Improperly labelled gene tree - no 'Biont' taxon.")

    # Improved version to allow different highest taxa, as long as they're identical:
    highest_taxa = [l.lineage[-1] for l in T]
    if len(set(highest_taxa)) != 1:
        raise ValueError("The highest taxa are not all equal. Please add the LCA taxon (e.g. 'biont') to the lineages.")
    return T


# Reconciliation functions:
def check_division(node, rank):
    """Checks if given ranks divides the children of the given annotated tree node
    in such a way that their evolutionary lineages are independent
    (i.e. there are no duplicate taxa below node).
    The lineages may contain 'UNKNOWN' value, which is ignored.

    Parameters
    ----------
    node : ete2.Tree object
        Lineage-annotated gene tree containing node to be checked. Annotation means that each leaf
        has to contain a lineage component that is a list of taxa representing
        the evolutionary lineage of the leaf. Lists have to be of equal length.
    rank : integer
        Index of the lineage list to be checked for proper lineage division.

    Value
    ----------
        True if division is proper, False if not.
    """
    assert isinstance(rank, int) & rank >= 0, "Rank should be a non-negative integer."
    if rank == 0:
        return True
    elif not node.children:
        return True
    else:
        tax_sets = [set(l.lineage[rank-1] for l in c) for c in node.children]
        tax_intersection = reduce(set.intersection, tax_sets)
        return not tax_intersection


def check_taxon_identity(node, rank):
    """Checks if species below node belong to the same taxa on the given taxonomic rank.

    Parameters
    ----------
    node : ete2.Tree object
        Lineage-annotated gene tree containing node to be checked. Annotation means that each leaf
        has to contain a lineage component that is a list of taxa representing
        the evolutionary lineage of the leaf. Lists have to be of equal length.
    rank : integer
        Index of the lineage list to be checked for taxon identity.

    Returns
    ----------
    out : bool
        True if leaf.lineage[rank] is equal for all leaves below node, False otherwise."""
    assert isinstance(rank, int) & rank >= 0, "Rank should be a non-negative integer."
    taxa = set(l.lineage[rank] for l in node)
    return len(taxa) == 1


def P(node, nb_of_ranks):
    """
    Computes the P characteristic of a node, i.e. the maximum taxonomic rank that satisfies the taxonomic separation
    of lineages.
    :param node: ete2.Tree
        A node of a gene tree.
    :param nb_of_ranks: int
        Length of a lineage.
    :return: int
        An index of a rank of the P characteristic.
    """
    for i in range(nb_of_ranks-1, -1, -1):
        divides = check_division(node, i)
        if divides:
            return i  # Note that rank 0 always separates.
    raise RuntimeError("No rank satisfies the lineage separation property. "
                       "This is impossible. Something went very wrong.")


def I(node, nb_of_ranks):
    """
    Computes the I characteristic of a node, i.e. the minimum taxonomic rank that satisfies the
    taxonomic identity of lineages.
    :param node: ete2.Tree
        Node of a gene tree.
    :param nb_of_ranks:
        Length of a lineage.
    :return: int
        An index of a rank of the I characteristic.
    """
    for i in range(nb_of_ranks):
        identity = check_taxon_identity(node, i)
        if identity:
            return i
    raise RuntimeError("No rank satisfies the identity. "
                       "Make sure that all the lineages have a common ancestor.")


def root_tree(tree, ranks):
    """Roots an unrooted, annotated tree based on evolutionary lineages."""
    tree = tree.copy()
    number_of_nodes = 0
    for node in tree.traverse():
        assert len(node.children) <= 3, "Input tree is not binary"  # i.e. no true multifurcations in gene tree
        number_of_nodes += 1
    rank_nmb = len(ranks)
    tree.unroot()
    old_outgroup = tree.children[0]
    outgroup_found = False
    i = 0
    while i <= number_of_nodes + 1:
        i += 1
        division_levels = [0 for _ in range(len(tree.children))]  # first levels that wrongly divide lineages when i-th branch is removed
        for c_id in range(len(tree.children)):
            binarized = tree.copy()
            binarized.remove_child(binarized.children[c_id])
            division = True
            r = -1
            while division and r < rank_nmb:
                r += 1
                division = check_division(binarized, r)
            division_levels[c_id] = r
        # new outgroup is chosen as the branch which after removal gives the lowest level
        # at which the lineages of the two remaining branches stop to be independent,
        # i.e. the lowest taxon that would be assigned to the top node of the remaining tree
        new_outgroup = tree.children[division_levels.index(min(division_levels))]
        # tree.set_outgroup(new_outgroup)
        # Check if we've returned to the same node as before:
        if new_outgroup == old_outgroup or new_outgroup.is_leaf():
            outgroup_found = True
            break
        # If not, the current root node becomes old_outgroup and then we set new_outgroup as the new root
        else:
            tree.remove_child(new_outgroup)
            old_outgroup = tree
            tree = new_outgroup
            tree.add_child(old_outgroup)
    tree.set_outgroup(old_outgroup)
    new_number_of_nodes = 0
    for node in tree.traverse():
        if len(node.children) > 2:
            raise RuntimeError("Error during evolutionary rooting: multifurcations introduced.\nOutput tree:\n%s", tree.write())
        new_number_of_nodes += 1
    if new_number_of_nodes > number_of_nodes:
        raise RuntimeError("Error during evolutionary rooting: new nodes introduced.\nOutput tree:\n%s", tree.write())
    elif new_number_of_nodes < number_of_nodes:
        raise RuntimeError("Error during evolutionary rooting: nodes lost.\nOutput tree:\n%s", tree.write())
    if not outgroup_found:
        raise ValueError("Could not find outgroup.\nOutput tree:\n%s" % tree.write())
    else:
        return tree


def get_most_common_taxon(node, rank):
    """Returns most common taxon occuring on a given rank under the given node."""
    children_taxa = [leaf.lineage[rank] for leaf in node]
    children_taxa.sort()
    children_taxa = [(t, sum(1 for _ in l) if t not in {"UNKNOWN", '0'} else 0) for t, l in groupby(children_taxa)]
    most_common_taxon = max(children_taxa, key=itemgetter(1))
    return most_common_taxon[0]


def label_nodes(tree, nb_of_ranks):
    """
    Labels nodes with the I and P characteristics.
    :param tree:
        Gene tree, labelled with the lineages.
    :param nb_of_ranks: int
        Number of evolutionary ranks, i.e. length of a lineage.
    :return: ete2.Tree
        A tree with I and P attributes at each node.
    """
    T = tree.copy()
    for node in T.traverse():
        Pv = P(node, nb_of_ranks)
        Iv = I(node, nb_of_ranks)
        node.P = Pv
        node.I = Iv
    return T


def is_transfer(node):
    """
    Returns true if the node is a possible transfer node.
    Assumes that P of parent of a transfer node has to be greater than 1
    (nodes with P==1 are duplication nodes).
    :param node: ete2.TreeNode
        A subtree annotated with I and P characteristics.
    :return: bool
    """
    if node.is_root():
        return False
    # elif node.up.P == 1:
    #     return False
    else:
        return node.P == node.I and node.I == node.up.I and node.up.P < node.up.I


def is_duplication(node):
    """
    Returns true if the node is a possible duplication node.
    Assumes that non-leaf nodes with P==1 are duplication nodes.
    :param node: ete2.TreeNode
        A subtree annotated with I and P characteristics.
    :return: bool
    """
    if node.is_leaf():
        return False
    elif node.P == 1:
        return True
    else:
        return len({node.I}.intersection({c.I for c in node.children})) != 0 and node.P < node.I


def classify_nodes(tree, ranks):
    """
    Classifies nodes into three groups: speciation, duplication and transfer nodes.
    If a node belongs to a duplication or transfer class, it indicates that a respective evolutionary event
    can be used to explain the topological discordance below this node. In other words, the tree is labelled
    with every possible evolutionary scenario that explains the discordance. The procedure does not decide which
    scenario to choose to obtain a parsimonious solution.
    :param tree:
        Annotated gene tree.
    :return: tree
        A tree with "event" attribute added to each leaf.
    """
    T = tree.copy()
    for node in T.traverse():
        if node.is_leaf():
            node.event = "Leaf"
            node.taxid = node.name.split('_')[0]
            node.rank = ranks[0]
            continue
        dpl = is_duplication(node)
        hgt = is_transfer(node)
        if dpl and hgt:
            node.event = "Duplication and Transfer"
            node.taxid = get_most_common_taxon(node, node.I)
            node.rank = ranks[node.I]
            raise RuntimeError("Duplication and Transfer occured. Check what happened.")
        elif dpl:
            node.event = "Duplication"
            node.taxid = get_most_common_taxon(node, node.I)
            node.rank = ranks[node.I]
        elif hgt:
            node.event = "Transfer"
            node.taxid = get_most_common_taxon(node, node.P)
            node.rank = ranks[node.P]
        else:
            node.event = "Speciation"
            assert node.I == node.P, "Improper speciation: P != I."
            node.taxid = get_most_common_taxon(node, node.I)
            node.rank = ranks[node.I]

    return T


def name_nodes(tree, taxon_to_name_mapping):
    """
    An attribute taxon_name is added to each node that has a taxid attribute.
    :param tree:
        Gene tree.
    :param taxon_to_name_mapping: dict
        A dictionary mapping taxIDs to taxon names.
    :return:
    """
    T = tree.copy()
    for node in T.traverse():
        try:
            node.scientific_name = taxon_to_name_mapping[node.taxid]
        except KeyError:
            node.scientific_name = 'UNKNOWN'
    return T


def process_tree(gene_tree, lineage_dict, ranks, rooted=False, max_unknown=3):
    """
    Labels, roots and reconciles a gene tree.
    :param gene_tree:
    :param lineage_dict:
    :return:
    """
    nb_of_ranks = len(ranks)
    assert nb_of_ranks > 0, "No ranks supplied!"
    check_data(lineage_dict, max_unknown)
    lineage_dict = curate_lineages(lineage_dict)
    gene_tree = annotate_tree(gene_tree, lineage_dict)
    check_annotation(gene_tree)
    if not rooted and len(gene_tree) > 2:
        midpoint = gene_tree.get_midpoint_outgroup()
        if midpoint:
            gene_tree.set_outgroup(midpoint)
        else:
            warn("No midpoint outgroup found! Attempting to proceed")
        try:
            gene_tree = root_tree(gene_tree, ranks=ranks)
        except:
            print "Error occured during rooting of the gene tree."
            raise
    gene_tree = label_nodes(gene_tree, nb_of_ranks)
    gene_tree = classify_nodes(gene_tree, ranks)
    return gene_tree


def is_acceptor(node):
    """
    Return true if the node is an acceptor node (i.e. edge from node.up to node is a transfer edge).
    Raises an error if the parent node is not a transfer node.
    :param node:
    :return: bool
    """
    if node.up.event != "Transfer":
        raise ValueError("Parent node is not a transfer node.")
    sibling = node.get_sisters()
    if len(sibling) != 1:
        raise RuntimeError("Improper number of sibling nodes")
    else:
        sibling = sibling[0]
    transfer_node = node.up
    transfer_sibling = transfer_node.get_sisters()
    if len(transfer_sibling) != 1:
        raise RuntimeError("Improper number of siblings of transfer node")
    else:
        transfer_sibling = transfer_sibling[0]
    node_I = -1
    node_tax_idty = False
    while not node_tax_idty:
        node_I += 1
        taxa = set([l.lineage[node_I] for l in node]).union(set([l.lineage[node_I] for l in transfer_sibling]))
        if len(taxa) == 1:
            node_tax_idty = True
    sibling_I = -1
    node_tax_idty = False
    while not node_tax_idty:
        sibling_I += 1
        taxa = set([l.lineage[sibling_I] for l in sibling]).union(set([l.lineage[sibling_I] for l in transfer_sibling]))
        if len(taxa) == 1:
            node_tax_idty = True
    if node_I > sibling_I:
        return True
    else:
        return False


def return_transfer_nodes(tree):
    """
    Returns a list of transfer nodes together with additional data.
    :param tree: ete2.Tree
        A processed gene tree.
    :return: list
        A list of tuples containing: transfer node, acceptor node, donor node and number of duplications required to explain
        the transfer in the DL model.
        Donor and acceptor nodes are nodes such that the transfer occured between edges directly above them.
    """
    transfer_nodes = [n for n in tree.traverse() if n.event == "Transfer"]
    data = []
    for n in transfer_nodes:
        # note that the root can't be a transfer node
        # leafs can't be transfer nodes as well because parents of leafs are always proper speciation nodes
        score = score_node(n)
        acceptor, donor = (n.children[0], n.children[1]) if is_acceptor(n.children[0]) else (n.children[1], n.children[0])
        data.append([n, acceptor, donor, score])
    return data


def name_internal_nodes(tree):
    tree = tree.copy()
    for node in tree.traverse():
        node.name = ''.join(n.name for n in node)
    return tree


def convert_to_lineages(tree):
    tree = tree.copy()
    _, nb_of_ranks = tree.get_farthest_leaf()
    nb_of_ranks = int(nb_of_ranks) + 1
    lineage_dict = dict()
    for leaf in tree:
        lineage = ['']*nb_of_ranks
        parent = leaf
        lineage[0] = leaf.name
        while not parent.is_root():
            parent = parent.up
            _, rank = parent.get_farthest_leaf()
            rank = int(rank)
            lineage[rank] = parent.name
        current_lineage = lineage[0]
        last_non_empty = 0
        for r in range(1, nb_of_ranks):
            if lineage[r] == '':
                lineage[r] = current_lineage + '*'*(r-last_non_empty)
            else:
                current_lineage = lineage[r]
                last_non_empty = r
        lineage_dict[leaf.name] = lineage
    return lineage_dict


def partition_tree(gene_tree, species_tree):
    """
    Partitions the gene tree, i.e. extracts independent evolutionary lineages.
    Returns a tree in which each node contains a mapping to a node in S, or None if the node has been removed.
    The mapping is a pointer to a node in species_tree.
    If the root of the species tree is not named, then the internal nodes are assumed to be not named as well.
    In this case, all internal nodes are automatically named by concatenating
    names of their leafs.
    :param gene_tree: ete2.Tree
        Tree to be partitioned. The names of the leafs have to correspond to the names of leafs of
        species_tree, possibly with identifier after an underscore (e.g. "a" or "a_1").
    :param species_tree: ete2.Tree
        Reference species tree.
    :return: tuple
        Gene tree with embedding into species tree and the resulting forest
    """
    gene_tree = gene_tree.copy()
    for l in gene_tree:
        l.name = l.name.split('_')[0]
    if not species_tree.name:
        species_tree = name_internal_nodes(species_tree)
    lineages = convert_to_lineages(species_tree)
    _, nb_of_ranks = species_tree.get_farthest_leaf()
    nb_of_ranks = int(nb_of_ranks) + 1
    leaf_mapping = dict([(l.name, species_tree & l.name) for l in gene_tree])
    gene_tree = annotate_tree(gene_tree, lineages)
    gene_tree = label_nodes(gene_tree, nb_of_ranks)
    for g in gene_tree.traverse():
        g.new_locus = False
    # list stored to obtain a mapping between tree and it's copy
    gene_tree_nodelist = [node for node in gene_tree.traverse(strategy="postorder")]
    working_tree = gene_tree.copy(method="deepcopy")
    for i, n in enumerate(working_tree.traverse(strategy="postorder")):
        n.origin_node = i
    not_embedded = [node for node in working_tree.traverse(strategy="postorder")]
    forest = []
    while not_embedded:
        working_tree.prune(not_embedded)
        working_tree = label_nodes(working_tree, nb_of_ranks)
        subtree = None
        for node in working_tree.traverse(strategy="postorder"):
            if not node.is_leaf() and node.I == 0:
                subtree = node
                break
            if not node.is_leaf() and node.I != node.P:
                subtree = node
                break
        if subtree and subtree.I == 0:  # profiting from lazy evaluation
            chosen_one = subtree.children[0]
        elif subtree:
            subtree_size = len(subtree)
            candidates = [n for n in subtree.children if n.I == node.I]
            candidates.extend([c for n in candidates for c in n.children])
            candidate_badness = [-1]*len(candidates)
            candidate_sizes = [len(c) for c in candidates]
            for i, c in enumerate(candidates):
                remaining_size = subtree_size - candidate_sizes[i]
                c_leafs = list(c)
                remaining_leafs = [l for l in subtree if l not in c_leafs]
                r_M = leaf_mapping[remaining_leafs[0].name].get_common_ancestor([leaf_mapping[l.name] for l in remaining_leafs])
                _, r_M_rank = r_M.get_farthest_leaf()
                r_M_rank = int(r_M_rank)
                # Check if the candidate can be removed, i.e. if no lineages mix after removal
                # Children of improper node can always be removed, but their children may not
                if c not in subtree.children:
                    leaf_set_1 = [l for l in remaining_leafs if l in subtree.children[0].get_leaves()]
                    leaf_set_2 = [l for l in remaining_leafs if l not in leaf_set_1]
                    tax_set_1 = set([l.lineage[r_M_rank - 1] for l in leaf_set_1])
                    tax_set_2 = set([l.lineage[r_M_rank - 1] for l in leaf_set_2])
                    if tax_set_1.intersection(tax_set_2):
                        # if the taxa intersect, the candidate can't be removed
                        continue
                c_M = leaf_mapping[c_leafs[0].name].get_common_ancestor([leaf_mapping[l.name] for l in c_leafs])
                candidate_diff = len(c_M) - candidate_sizes[i]
                second_subtree_diff = len(r_M) - remaining_size
                assert candidate_diff >= 0 and second_subtree_diff >= 0, "Negative score! Investigate"
                candidate_badness[i] = candidate_diff + second_subtree_diff
            max_badness = max(candidate_badness)
            # candidate badness is non-negative, so -1 means that node can't be removed
            candidate_badness = [c if c != -1 else max_badness + 9000 for c in candidate_badness]  # IT'S OVER 9000!
            chosen_one = candidates[candidate_badness.index(min(candidate_badness))]
        else:
            chosen_one = working_tree
        gene_tree_nodelist[chosen_one.origin_node].new_locus = True
        for n in chosen_one.traverse():
            n_leafs = list(n)
            gene_tree_nodelist[n.origin_node].add_feature("embedding", leaf_mapping[n_leafs[0].name].get_common_ancestor([leaf_mapping[l.name] for l in n_leafs]))
        forest.append(chosen_one)
        working_tree = gene_tree.copy(method="deepcopy")
        for i, n in enumerate(working_tree.traverse(strategy="postorder")):
            n.origin_node = i
        not_embedded = [node for node in working_tree if not hasattr(node, "embedding")]
    return (gene_tree, forest)


def transfer_score(node):
    """
    Returns estimated number of duplications needed to explain the incongruence with the DL model.
    If the LCA node of transfer donor and acceptor is not in the gene tree, 9001 is returned.
    :param node
    :return: int
    """
    node_taxa = set([l.name.split('_')[0] for l in node])
    parent = node
    score = 0
    while not parent.is_root():
        parent = parent.up
        score += 1
        cdrn = parent.children
        cdrn_txs = [set([l.name.split('_')[0] for l in c]) for c in cdrn]
        intrsct = [node_taxa.intersection(c_tx) for c_tx in cdrn_txs]
        if all(intrsct):
            return score
    return 9001


if __name__=="__main__":
    S = ete2.Tree('(((a, b), (c, d)), e, (f, g));')
    G = ete2.Tree('((((a, c), b), (e, f)), (c, ((g, f), (d, a))));')
    S = ete2.Tree('(a, b);')
    G = ete2.Tree('((a, b), (a, b));')
    # S = ete2.Tree('((v, s), (x, w));')
    # G = ete2.Tree('(((s, v), x), (w, w));')
    # universe = list(string.ascii_lowercase)
    # S = ete2.Tree()
    # S.populate(20, names_library=universe)
    # species_names = [l.name for l in S]
    # G = ete2.Tree()
    # G.populate(30, names_library=species_names, reuse_names=True)
    print S
    print G

    # Not needed anymore as performed in partition_tree:
    # S = name_internal_nodes(S)
    # lgs = convert_to_lineages(S)
    # _, nb_of_ranks = S.get_farthest_leaf()
    # nb_of_ranks = int(nb_of_ranks) + 1
    # G = annotate_tree(G, lgs)
    # G = label_nodes(G, nb_of_ranks)
    g, forest = partition_tree(G, S)
    for n in g.traverse():
        if hasattr(n, 'embedding'):
            n.embed_name = n.embedding.name
        else:
            n.embed_name = 'N/A'
        if n.new_locus:
            n.score = transfer_score(n)
        else:
            n.score = 0
    print g.get_ascii(show_internal=True, attributes=['embed_name', "new_locus", "score"])
    for f in forest:
        print f.get_ascii(show_internal=True)

    # Porownanie z algorytmem optymalnym
    # reps = 100
    # gene_tree_sizes = list(range(21))
    # # forest_sizes = [0]*len(gene_tree_sizes)
    # forest_sizes = np.array([[0] * (reps+1)] * len(gene_tree_sizes))
    # forest_sizes[:, 0] = gene_tree_sizes
    # start_time = time()
    # for i, size in enumerate(gene_tree_sizes):
    #     # tmp_sizes = [0]*reps
    #     print "Gene tree size", size
    #     for j in range(reps):
    #         S = ete2.Tree()
    #         S.populate(10, names_library=universe)
    #         species_names = [l.name for l in S]
    #         G = ete2.Tree()
    #         G.populate(size, names_library=species_names, reuse_names=True)
    #         g, forest = partition_tree(G, S)
    #         forest_sizes[i, j+1] = len(forest)
    #     #     tmp_sizes[j] = len(forest)
    #     # forest_sizes[i] = sum(tmp_sizes)/float(reps)
    # end_time = time()
    # runtime = end_time - start_time
    #
    # np.savetxt("Heur_forest_sizes.tsv", forest_sizes, delimiter='\t', fmt='%i')
    #
    # average_sizes = np.mean(forest_sizes[:,1:], axis = 1)
    # plt.plot(gene_tree_sizes, average_sizes)
    # plt.title("Average forest size from %i replications.\nHeuristic algorithm. Runtime: %.2f" % (reps, runtime))
    # plt.xlabel("Number of leafs in G")
    # plt.ylabel("Forest size")
    # plt.savefig("Heuristic_average_size.png")
    # plt.show()

    # Ekstremalne rozmiary drzew
    # reps = 10
    # gene_tree_sizes = list(range(0, 601, 100))
    # # forest_sizes = [0]*len(gene_tree_sizes)
    # forest_sizes = np.array([[0] * (reps+1)] * len(gene_tree_sizes))
    # forest_sizes[:, 0] = gene_tree_sizes
    # start_time = time()
    # for i, size in enumerate(gene_tree_sizes):
    #     # tmp_sizes = [0]*reps
    #     print "Gene tree size", size
    #     for j in range(reps):
    #         S = ete2.Tree()
    #         S.populate(200)
    #         species_names = [l.name for l in S]
    #         G = ete2.Tree()
    #         G.populate(size, names_library=species_names, reuse_names=True)
    #         g, forest = partition_tree(G, S)
    #         forest_sizes[i, j+1] = len(forest)
    #     #     tmp_sizes[j] = len(forest)
    #     # forest_sizes[i] = sum(tmp_sizes)/float(reps)
    # end_time = time()
    # runtime = end_time - start_time
    #
    # np.savetxt("Extreme_forest_sizes.tsv", forest_sizes, delimiter='\t', fmt='%i')
    #
    # average_sizes = np.mean(forest_sizes[:,1:], axis = 1)
    # plt.plot(gene_tree_sizes, average_sizes)
    # plt.title("Average forest size from %i replications.\nHeuristic algorithm. Runtime: %.2f" % (reps, runtime))
    # plt.xlabel("Number of leafs in G")
    # plt.ylabel("Forest size")
    # plt.savefig("Extreme_average_size.png")
    # plt.show()

















