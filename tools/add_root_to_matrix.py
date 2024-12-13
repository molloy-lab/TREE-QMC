import argparse
import numpy
import sys


def add_root_to_matrix(fin, fout):
    # Read first row with mutation names
    line = fin.readline()
    fout.write(line)

    data = line.split()
    muts = data[1:]
    nmut = len(muts)

    for line in fin:
        fout.write(line)

    fout.write("root\t")
    fout.write('\t'.join(['0'] * nmut))


def main(args):
    with open(args.input, 'r') as fin, \
         open(args.output, 'w') as fout:
        add_root_to_matrix(fin, fout)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input file (with matrix)",
                        required=True)

    parser.add_argument("-o", "--output", type=str,
                        help="Output prefix",
                        default="")

    main(parser.parse_args())
