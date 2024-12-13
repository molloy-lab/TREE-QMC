import argparse
import numpy
import treeswift
import sys


def matrix_to_trees(fin):
    # Read header for number of mutations
    aboves = []
    belows = []
    line = fin.readline()
    muts = line.split()
    for mut in muts[1:]:
        aboves.append([])
        belows.append([])
    nmut = len(aboves)
        
    # Read matrix data
    for line in fin:
        data = line.split()
        cell = data[0]
        for i, x in enumerate(data[1:]):
            if x == '0':
                aboves[i].append(cell)
            elif x == '1':
                belows[i].append(cell)

    # Create trees (newick strings)
    trees = []
    for i in range(nmut):
        above = aboves[i]
        below = belows[i]

        na = len(above)
        nb = len(below)

        if na + nb > 3:
            if nb > 1:
                astr = ','.join(above)
                bstr = '(' + ','.join(below) + ')'
                tree = '(' + astr + ',' + bstr + ");\n"
                trees.append(tree)

    return trees


def main(args):
    with open(args.input, 'r') as fin:
        trees = matrix_to_trees(fin)

    with open(args.output, 'w') as fout:
        for tree in trees:
            fout.write(tree)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input file (with matrix)",
                        required=True)

    parser.add_argument("-o", "--output", type=str,
                        help="Output prefix",
                        default="")

    main(parser.parse_args())
