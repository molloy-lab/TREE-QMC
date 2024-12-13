import argparse
import sys


def huntress_matrix_to_scistreemat(fin, fout, a, b):
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
        cdata.append(data[1:])

    n = len(cnams)
    m = len(mnams)

    fout.write("HAPLOID %d %d " % (m, n))
    fout.write(' '.join(cnams))
    fout.write('\n')
    for i in range(m):
        fout.write(mnams[i])
        for j in range(n):
            x = cdata[j][i]
            if x == '0':
                fout.write(' ' + b)
            elif x == '1':
                fout.write(' ' + a)
            elif x == '3':
                fout.write(' 0.5')
            else:
                sys.exit("Unknown matrix element!")
        fout.write('\n')


def main(args):
    a = str("%1.12f" % args.alpha)
    b = str("%1.12f" % (1 - args.beta))

    with open(args.input, 'r') as fin, \
         open(args.output, 'w') as fout:
        huntress_matrix_to_scistreemat(fin, fout, a, b)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input file (with matrix)",
                        required=True)

    parser.add_argument("-a", "--alpha", type=float,
                        help="False positive rate",
                        required=True)

    parser.add_argument("-b", "--beta", type=float,
                        help="False negative rate",
                        required=True)

    parser.add_argument("-o", "--output", type=str,
                        help="Output prefix",
                        default="")

    main(parser.parse_args())
