from lineage_separation import partition_tree, transfer_score, name_internal_nodes, convert_to_lineages, annotate_tree, label_nodes
import getopt
import sys
import ete2

help_text = """
NAME
    decompose_history.py
USAGE
    python decompose_history.py -g GENE_TREE -s SPECIES_TREE -f
DESCRIPTION
    Decomposes the evolutionary history of loci in the gene tree, based on
    the supplied species tree.
    If the species tree's internal nodes are named, then those names are used.
    Otherwise, new names are created by concatenating names of leafs.
OPTIONS
    -f: Return forest
        Returns the forest obtained after decomposition.
        Otherwise, returns a gene tree with NHX attributes containing indicator if a node
        is the first observation of a new loci (NHX new_locus), the embedding
        of the node into the species tree given by the species tree's node name (NHX embedding),
        and the post-hoc transfer score (NHX score).

"""


if __name__ == "__main__":
    return_forest = False
    species_tree = None
    gene_tree = None
    opts, args = getopt.getopt(sys.argv[1:], "g:s:f")
    for opt, arg in opts:
        if opt == '-g':
            gene_tree = arg
        elif opt == '-s':
            species_tree = arg
        elif opt == '-f':
            return_forest = True

    if gene_tree[-1] != ';':
        gene_tree += ';'
    if species_tree[-1] != ';':
        species_tree += ';'

    species_tree = ete2.Tree(species_tree, format=8)
    gene_tree = ete2.Tree(gene_tree, format=2)

    gene_tree, forest = partition_tree(gene_tree, species_tree)

    for n in gene_tree.traverse():
        if hasattr(n, 'embedding'):
            n.embedding = n.embedding.name
        else:
            n.embedding = 'N/A'
        if n.new_locus:
            n.score = transfer_score(n)
        else:
            n.score = 0

    if return_forest:
        print ' '.join(f.write(format=9) for f in forest)
    else:
        print gene_tree.write(format=2, features=["embedding", "new_locus", "score", "I", "P"], format_root_node=True)

