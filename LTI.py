#! /usr/bin/python
from src import heuristic, dynamic, multitree
from src.generic import Decomposition
import getopt
import sys
import ete3

help_text = """
NAME
    LTI: Locus Tree Inference
USAGE
    LTI [OPTIONS] 
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
    -s: Species tree
    -g: Gene tree
        Either paths to files with the trees in Newick format, or the trees itself as strings.
        The tree in Newick format needs to end with a semicolon.
        NHX annotations are allowed, but will not be preserved in the output.
        The gene tree leaves need to follow a SPECIES_GENE convention,
        i.e. each leaf name of the gene tree needs to contain the name of the corresponding
        species, followed by an underscore, followed by an identifier.
    -a: Algorithm used
        The algorithm used for decomposition. Possible options:
            D: dynamic programming; optimal decomposition in quadratic time and memory
            H: heuristic; less optimal results in linear memory 
            A: automatic choice; dynamic programming for gene trees with less than 40 leafs.
        Default is the automatic choice. 
    -f: Return forest
        Returns the forest obtained after decomposition.
        Otherwise, returns a gene tree with NHX attributes containing indicator if a node
        is the first observation of a new loci (NHX new_locus), the embedding
        of the node into the species tree given by the species tree's node name (NHX embedding),
        and the post-hoc transfer score (NHX score).
    -r: Return roots
        Returns the roots file. This is a file with three lines. The first one contains the gene tree,
        the second one the species tree, and the third one the list of postorder IDs of root nodes 
        from the decomposition forest.
    -e: Gene tree format
        Newick format of the input gene tree (as in ete2). Default: 9, leaf names only.
        The format will be preserved in the output. 
    -n: Suppress automatic ranking
        If this flag is not set, the species tree will be ranked automatically based on the nodes' topology.
        Alternatively, one may specify the ranks as NHX attributes. Example: 
        '(S3[&&NHX:rank=0],(S2[&&NHX:rank=0],S1[&&NHX:rank=0])[&&NHX:rank=1]);'
        Currently, custom ranks are supported only in heuristic algorithm.
    
"""


if __name__ == "__main__":
    return_format = "tree"
    species_tree = None
    gene_tree = None
    gene_format = 9
    automatic_ranking = True
    dynamic_programming = True
    select_algorithm = True

    if len(sys.argv) <= 1:
        print(help_text)
        quit()
    opts, args = getopt.getopt(sys.argv[1:], "g:s:e:fhrna:")
    for opt, arg in opts:
        if opt == '-g':
            gene_tree = arg
        elif opt == '-s':
            species_tree = arg
        elif opt == '-f':
            return_format = "forest"
        elif opt == '-r':
            return_format = "roots"
        elif opt == '-h':
            print(help_text)
            quit()
        elif opt == '-e':
            gene_format = int(arg)
        elif opt == '-n':
            automatic_ranking = False
        elif opt == '-a':
            if arg.upper() == 'H':
                dynamic_programming = False
                select_algorithm = False
            elif arg.upper() == "D":
                select_algorithm = False
            elif arg.upper() != "A":
                raise ValueError("Selected algorithm not recognized.")

    # using ete's parser to get info on trees
    S = ete3.Tree(species_tree, format=8)
    G = ete3.Tree(gene_tree, format=gene_format)
    species_names = set(l.name for l in S)
    original_gene_names = [l.name for l in G]
    for i, l in enumerate(G):
        l.name = '_'.join(l.name.split('_')[:-1])
        assert (l.name in species_names), 'Gene tree leaf %s does not correspond to any species!' % original_gene_names[i]
    
    if automatic_ranking:
        from src.generic import assign_ranks
        assign_ranks(S)
        
    # Decomposing
    if not dynamic_programming or len(G) > 20:
        roots = heuristic.decompose(G, S)
    else:
        Gc = G.copy('deepcopy')
        MS = S.write(format=9)
        MG = Gc.write(format=9)
        MS = MS.replace(';', '')
        MG = MG.replace(';', '')
        MS = multitree.Tree(MS)
        MG = multitree.Tree(MG)
        roots = dynamic.decompose(MG, MS)
        roots = next(roots)
    # Restoring original leaf names:
    for on, l in zip(original_gene_names, G):
        l.name = on

    # Returning results
    D = Decomposition(G, S, roots)
    if return_format == "forest":
        forest = D.forest()
        print(' '.join(f.write(format=gene_format) for f in forest))
    elif return_format == "tree":
        locus_tree = D.locus_tree()
        D._assert_sources()
        print(locus_tree.write(format=gene_format, features=["source", "locus_id", "I", "P"], format_root_node=True))
    elif return_format == "roots":
        for g in G:
            g.name = g.name.split('_')[0]
        to_write = G.write(format=gene_format).replace(';', '') + '\n'
        to_write += S.write(format=9).replace(';', '') + '\n'
        to_write += '(' + ','.join(map(str, D.roots)) + ')'
        print(to_write)

