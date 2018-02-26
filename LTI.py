from src.lineage_separation import partition_tree, transfer_score
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
    Species tree needs to contain only the names of nodes (no supports, branches etc.)
    By default, this requirements holds also for gene trees. If you wish to
    preserve branch lenghts and/or supports in the results, please
    indicate the input gene tree format (as specified by ete3 documentation). 
OPTIONS
    -f: Return forest
        Returns the forest obtained after decomposition.
        Otherwise, returns a gene tree with NHX attributes containing indicator if a node
        is the first observation of a new loci (NHX new_locus), the embedding
        of the node into the species tree given by the species tree's node name (NHX embedding),
        and the post-hoc transfer score (NHX score).
    -e: Gene tree format
        Newick format of the input gene tree. Default: 9 (leaf names only)

"""


if __name__ == "__main__":
    return_forest = False
    species_tree = None
    gene_tree = None
    gene_format = 9
    if len(sys.argv) <= 1:
        print help_text
        quit()
    opts, args = getopt.getopt(sys.argv[1:], "g:s:e:fh")
    for opt, arg in opts:
        if opt == '-g':
            gene_tree = arg
        elif opt == '-s':
            species_tree = arg
        elif opt == '-f':
            return_forest = True
        elif opt == '-h':
            print help_text
            quit()
        elif opt == '-e':
            gene_format = int(arg)

    # if gene_tree[-1] != ';':
    #     gene_tree += ';'
    # if species_tree[-1] != ';':
    #     species_tree += ';'

    species_tree = ete2.Tree(species_tree, format=8)
    gene_tree = ete2.Tree(gene_tree, format=gene_format)
    
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
        print gene_tree.write(format=gene_format, features=["embedding", "new_locus", "score", "I", "P"], format_root_node=True)

