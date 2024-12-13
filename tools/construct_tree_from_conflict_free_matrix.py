import argparse
import numpy
import treeswift
import sys


def build_tree_from_clades(taxa, clades):
    # Get unique clades
    clade_list = {}
    for clade in clades:
        if len(clade) > 0:
            cid = ','.join(sorted(clade))
            clade_list[cid] = clade

    # Build star tree
    tree = treeswift.Tree()
    root = None
    for node in tree.traverse_preorder():
        if node.is_root():
            root = node

    tax2nod = {}
    for x in taxa:
        node = treeswift.Node()
        node.set_label(x)
        root.add_child(node)
        tax2nod[x] = node

    # Sort clades by size
    cids = [cid for cid in clade_list.keys()]
    sizes = []
    for cid in cids:
        sizes.append(len(clade_list[cid]))

    sizes = numpy.array(sizes)
    inds = numpy.argsort(sizes)
    inds = numpy.flip(inds)

    # Refine clades largest to smallest
    for i in inds:
        clade = clade_list[cids[i]]

        leaf = tax2nod[clade[0]]
        old_parent = leaf.get_parent()
        new_parent = treeswift.Node()
        old_parent.add_child(new_parent)

        for tax in clade:
            leaf = tax2nod[tax]
            parent = leaf.get_parent()

            if parent != old_parent:
                sys.exit("Not conflict free!")

            old_parent.remove_child(leaf)
            new_parent.add_child(leaf)

    return tree


def matrix_to_clades(fin):
    clades = []
    
    # Read header
    line = fin.readline()

    muts = line.split()
    for mut in muts[1:]:
        clades.append([])

    taxon_set = set([])
        
    # Read matrix
    for line in fin:
        data = line.split()
        cell = data[0]
        for i, x in enumerate(data[1:]):
            if float(x) == 1.0:
                clades[i].append(cell)
            taxon_set.add(cell)

    taxa = list(taxon_set)
    return [taxa, clades]


def main(args):
    with open(args.input, 'r') as fin:
        [taxa, clades] = matrix_to_clades(fin)

    tree = build_tree_from_clades(taxa, clades)

    with open(args.output, 'w') as fout:
        fout.write(tree.newick())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input file",
                        required=True)

    parser.add_argument("-o", "--output", type=str,
                        help="Output prefix",
                        default="")

    main(parser.parse_args())
