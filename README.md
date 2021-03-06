# Locus Tree Inference
 
Locus Tree Inference in the parsimony framework is a method of decomposing a gene tree into a set of subtrees such that each subtree is embeddable into the species tree. 
Each subtree represents the evolutionary history of a single locus. 

The software can be used, among others, to:
* Identify sources of incongruences between gene and species trees;
* Indentify evolutionary events, like the Horizontal Gene Transfer or Gene Duplication;
* Quantify the level of incongruence between a pair of trees.

Currently, the gene tree needs to be binary, but the species tree can contain polytomies.
The Locus Tree Inference method uses only the information about the trees' topologies (i.e. no branch lengths or supports needed). 

Note that this is still a work in progress, so if you notice any inconvenience in using the software or want us to add some functionality, do not hesitate to leave a comment or contact us at m_ciach@student.uw.edu.pl.

The LTI software is developed using Python 3 and the [ETE toolkit](http://etetoolkit.org/).

## Installation

To start decomposing your trees, simply clone the repository by typing into the commandline:

```shell 
git clone https://github.com/mciach/LocusTreeInference
```

## Quick start 

The basic usage of the software is 

```shell
python3 LTI.py -f -g GENE TREE -s SPECIES TREE 
```

where `GENE TREE` and `SPECIES TREE` are strings encoding trees in Newick format. 
The names of leaves of the have to correspond to the names of leaves of the species tree, possibly with an identifier after an underscore.
Running this command will return a forest of subtrees resulting from decomposition of the GENE TREE.

### Examples

```shell 
python3 LTI.py -f -g '(((a, c), d), (a, b));' -s '((a, b), (c, d));'
```

```shell 
python3 LTI.py -f -g '(((a_1, c_1), d_1), (a_2, b_1));' -s '((a, b), (c, d));'
```

## Usage

The general syntax for the software is 

```shell 
python3 LTI.py [OPTIONS]
```

where the possible `OPTIONS` include:

* `-h`: print the help message and exit;
* `-g`: specify the gene tree;
* `-s`: specify the species tree; 
* `-f`: return the result as a forest;

If the flag `-f` is not specified, the locus tree is returned as a Newick string with NHX annotations. 
The annotation is designed for the visualization program (visualize_decomposition.py). 

Additional options are described in the help message.

## Licensing and Credit

If you use LTI in your work, please cite:

Ciach, Michał Aleksander, Anna Muszewska, and Paweł Górecki. "Locus-aware decomposition of gene trees with respect to polytomous species trees." Algorithms for Molecular Biology 13.1 (2018): 11.

