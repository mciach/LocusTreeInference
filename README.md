# LocusTreeInference
Locus Tree Inference in parsimony framework.

## USAGE
>    python decompose_history.py -g GENE_TREE -s SPECIES_TREE -f
## DESCRIPTION
    Decomposes the evolutionary history of loci in the gene tree, based on
    the supplied species tree (both in Newick format).
    If the species tree's internal nodes are named, then those names are used.
    Otherwise, new names are created by concatenating names of leafs.
    Currently, the ranks for the species tree are created automatically. 
## OPTIONS
    -f: Return forest
        Returns the forest obtained after decomposition.
        Otherwise, returns a gene tree with NHX attributes containing indicator if a node
        is the first observation of a new loci (NHX new_locus), the embedding
        of the node into the species tree given by the species tree's node name (NHX embedding),
        and the post-hoc transfer score (NHX score).
