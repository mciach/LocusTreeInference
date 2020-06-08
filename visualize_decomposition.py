import getopt
import sys
from ete3 import Tree, NodeStyle, faces, TreeStyle
from ete3 import NCBITaxa
# print Tableau_20.hex_colors
colors = ["#a0a0af", "#a8cfe6", "#237ab7", "#b3e08c", "#35aa2f", "#fe9d9c",
                  "#f72d2f", "#fec271", "#ff8010", "#cdb4d8"]

EVOLUTIONARY_RANKS = ['species', 'genus', 'family', 'order', 'class', 'phylum']

leaf_name_speparation = 1

doc = """NAME:
    visualize_decomposition.py 
USAGE:
    python visualize_decomposition.py [OPTIONS] FILE
DESCRIPTION:
    Visualize the decomposed gene tree. The decomposition needs to be provided in FILE,
    represented as NHX-annotated gene tree (i.e. the default output of LTI.py).
    The different loci are highlighted by different colours. Numbers above branches
    denote the values of I/P mappings, numbers below branches denote branch supports (if enabled).   
    The mapping values are parsed from FILE.
OPTIONS:
    -h
        Print this message and exit.
    -e: integer from 0 to 9, default 9
        The format of the input locus tree (see the documentation of ETE toolkit). 
        This is useful mostly in the case when the application cannot parse
        the supplied gene tree; in such case, adjusting the format usually solves the problem.
    -u: integer, 0 or 1, default 1
        Whether to make the locus tree ultrametric (1) or preserve original branch lengths (0). 
        Default is 1, but only because in our opinion it looks better. Feel free to disable it.
    -s: float in [0, 1], default 1
        Branch support threshold. Only supports above this value will be shown. 
        Set to 1 to disable showing supports, set to 0 to show all supports.
    -i: integer, 0 or 1, default 0
        If set to 1, node IDs in postorder numbering are rendered.
"""

if __name__=="__main__":
    # species_tree = None
    gene_tree = None
    gene_format = 9
    ultrametric = 1
    support_threshold = 1.
    show_IDs = 0
    additional_features = []
    opts, args = getopt.getopt(sys.argv[1:], "hg:e:a:u:s:i:")
    for opt, arg in opts:
        if opt == "-e":
            gene_format = int(arg)
        elif opt == "-a":
            additional_features = arg.split()
        elif opt == '-h' or not args:
            print(doc)
            quit()
        elif opt == '-u':
            ultrametric = int(arg)
        elif opt == '-s':
            support_threshold = float(arg)
        elif opt == '-i':
            show_IDs = int(arg)
    gene_tree = args[0]
    print(gene_tree)
    print(gene_format)
    gene_tree = Tree(gene_tree, format=gene_format)


    #gene_tree = ete2.Tree("/home/ciach/Projects/Tree_node_classification/Comparison_with_Notung/Sequences/ABV48733.1/decomposed_tree.dnd")
    #species_tree = ete2.Tree("/home/ciach/Projects/Tree_node_classification/Comparison_with_Notung/Sequences/ABV48733.1/species_tree.dnd")

    suppFace = faces.AttrFace("support", text_suffix=" ", formatter="%.2f", fsize=8)
    iFace = faces.AttrFace("I", fsize=8, text_suffix='/')
    pFace = faces.AttrFace("P", fsize=8, text_suffix=' ')
    idFace = faces.AttrFace("id", fsize=8)
    for i, g in enumerate(gene_tree.traverse(strategy="postorder")):
        g.id = i
        #g.name = g.name.split('_')[0]
        g.coloured = False
        g.nw = g.source
        g.source = True if g.source == 'True' else False
        g.dist = 1.0
        if hasattr(g, 'I'):
            if g.I == "None":
                g.I = -1
            else:
                g.I = int(g.I) + 1
        else:
            g.I = -1
        if hasattr(g, 'P'):
            if g.P == "None":
                g.P = -1
            else:
                g.P = int(g.P) + 1
        else:
            g.P = -1
        if not g.is_leaf():
            if float(g.support) > support_threshold:
                g.add_face(suppFace, position="branch-bottom", column=-1)
            g.add_face(iFace, position="branch-top", column=-1)
            g.add_face(pFace, position="branch-top", column=0)
        if show_IDs:
            g.add_face(idFace, position="branch-right", column=1)

    current_color = -1

    for g in gene_tree.traverse(strategy="postorder"):
        g.name = ' '*leaf_name_speparation + g.name + ' '*leaf_name_speparation
        if g.source:
            current_color += 1
            current_color %= len(colors)
            style = NodeStyle()
            style['vt_line_color'] = colors[current_color]
            style['hz_line_color'] = colors[current_color]
            style['size'] = 0
            style['fgcolor'] = '#000000'
            style["vt_line_width"] = 2
            style["hz_line_width"] = 2
            for gg in g.traverse():
                if not gg.coloured:
                    gg.set_style(style)
                    gg.coloured = True
            source_style = NodeStyle(style)
            source_style['size'] = 5
            source_style['fgcolor'] = '#C01000'
            g.set_style(source_style)

    tstyle = TreeStyle()
    #tstyle.show_leaf_name = False
    tstyle.scale = 28
    tstyle.branch_vertical_margin = 6
    tstyle.show_branch_length = False
    
    # tstyle.show_branch_support = True
    tstyle.show_scale = False

    if ultrametric == 1:
        gene_tree.convert_to_ultrametric()
    gene_tree.show(tree_style=tstyle)
    # species_tree.show()
