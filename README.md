# LocusTreeInference
Locus Tree Inference in parsimony framework.

## USAGE
>    python decompose_history.py -g GENE_TREE -s SPECIES_TREE -f
## DESCRIPTION
    Decomposes the evolutionary history of loci in the gene tree, based on
    the supplied species tree (both in Newick format).
    The names of the leaves of the gene tree have to correspond to the names of the leaves of the species tree.
## EXAMPLE
>   python decompose_history.py -f -g '(((a, c), d), (a, b));' -s '((a, b), (c, d));'
## OPTIONS
    -f: Return forest
        Returns the forest obtained after decomposition.
        Otherwise, returns a gene tree with NHX attributes containing indicator if a node
        is the first observation of a new loci (NHX new_locus), the embedding
        of the node into the species tree given by the species tree's node name (NHX embedding),
        and the post-hoc transfer score (NHX score).
