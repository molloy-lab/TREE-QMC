import argparse
import sys


def huntress_matrix_to_fastme(fin, fout):
    # Read first row with mutation names
    line = fin.readline()
    data = line.split()
    mnams = data[1:]
        
    # Read rest of rows
    cnams = []
    cdata = []
    for line in fin:
        data = line.split()
        cnams.append(data[0])

        cdat = []
        for x in data[1:]:
            if x == '0':
                cdat.append('A')
            elif x == '1':
                cdat.append('T')
            else:
                cdat.append('-')

        cdata.append(''.join(cdat))

    n = len(cnams)
    m = len(mnams)

    fout.write("%d %d\n" % (n, m))
    for cnam, cdat in zip(cnams, cdata):
        fout.write(cnam + ' ' + cdat + '\n')


def main(args):
    with open(args.input, 'r') as fin, \
         open(args.output, 'w') as fout:
        huntress_matrix_to_fastme(fin, fout)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input file (with matrix)",
                        required=True)

    parser.add_argument("-o", "--output", type=str,
                        help="Output prefix",
                        default="")

    main(parser.parse_args())
