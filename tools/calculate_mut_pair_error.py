import argparse
import numpy
import sys


def read_matrix(fin):
    # Read header for number of mutations
    line = fin.readline()
    data = line.split()

    muts = data[1:]
    nmut = len(muts)

    cols = []
    for i in range(nmut):
        cols.append([])

    # Read cell rows and transpose so mutations are rows
    cells = []
    for line in fin:
        data = line.split()
        cell = data[0]

        cells.append(cell)
        for i in range(nmut):
            cols[i].append(float(data[i+1]))

    mcell = len(cells)

    matrix = {}
    for i in range(nmut):
        matrix[muts[i]] = numpy.array(cols[i])

    return matrix


def get_mutation_relationship(mut_i, mut_j):
    total = len(mut_i)
    if total != len(mut_j):
        sys.exit("Mutations are different sizes!\n")

    i0_j0 = 0
    i0_j1 = 0
    i1_j0 = 0
    i1_j1 = 0
    for i, j in zip(mut_i, mut_j):
        if i == 0:
            if j == 0:
                i0_j0 += 1
            elif j == 1:
                i0_j1 += 1
        else:
            if j == 0:
                i1_j0 += 1
            elif j == 1:
                i1_j1 += 1

    same = i0_j0 + i1_j1

    if same == total:
        # Mutation i and j occurred on same branch
        return 0
    elif same + i0_j1 == total:
        # Mutation i occurred on branch ancestral of mutation j
        return 2
    elif same + i1_j0 == total:
        # Mutation j occurred on branch ancestral of mutation i
        return 3
    elif i0_j0 + i0_j1 + i1_j0 == total:
        # Mutation i,j occurred on different lineages
        return 1
    else:
        sys.exit("Mutations are incompatible!\n")


def update_error(error, true_ij, esti_ij):
    # If mutations are on the same lineage
    if true_ij == 0:
        if esti_ij == 0:
            error["SL"]["TP"] += 1
        else:
            error["SL"]["FN"] += 1
    else:
        if esti_ij == 0:
            error["SL"]["FP"] += 1
        else:
            error["SL"]["TN"] += 1


    # If mutations are on different lineages
    if true_ij == 1:
        if esti_ij == 1:
            error["DL"]["TP"] += 1
        else:
            error["DL"]["FN"] += 1
    else:
        if esti_ij == 1:
            error["DL"]["FP"] += 1
        else:
            error["DL"]["TN"] += 1

    # If mutations are ancestor-descendant
    if true_ij == 2:
        if esti_ij == 2:
            error["AD"]["TP"] += 1
        elif esti_ij == 3:
            error["AD"]["FLIP"] += 1
        else:
            error["AD"]["FN"] += 1
    elif true_ij == 3:
        if esti_ij == 2:
            error["AD"]["FLIP"] += 1
        elif esti_ij == 3:
            error["AD"]["TP"] += 1
        else:
            error["AD"]["FN"] += 1
    else:
        if esti_ij == 2:
            error["AD"]["FP"] += 1
        elif esti_ij == 3:
            error["AD"]["FP"] += 1
        else:
            error["AD"]["TN"] += 1


def calculate_error(true_mat, esti_mat):
    error = {}
    for x in ["SL", "DL", "AD"]:
        error[x] = {}
        error[x]["TP"] = 0
        error[x]["TN"] = 0
        error[x]["FP"] = 0
        error[x]["FN"] = 0
    error["AD"]["FLIP"] = 0

    true_muts = set([k for k in true_mat.keys()])
    esti_muts = set([k for k in esti_mat.keys()])

    muts = list(true_muts.intersection(esti_muts))
    nmut = len(muts)

    #total = 0
    for i in range(0, nmut - 1):
        mi = muts[i]
        true_mut_i = true_mat[mi]
        esti_mut_i = esti_mat[mi]
        for j in range(i + 1, nmut):
            mj = muts[j]
            true_mut_j = true_mat[mj]
            esti_mut_j = esti_mat[mj]

            true_ij = get_mutation_relationship(true_mut_i, true_mut_j)
            esti_ij = get_mutation_relationship(esti_mut_i, esti_mut_j)

            update_error(error, true_ij, esti_ij)
            #total += 1

    #real_total = (nmut * (nmut - 1)) / 2
    #if (total != real_total):
    #    sys.exit("Did not test all pairs!\n")

    line = ""
    for x in ["SL", "DL", "AD"]:
        tp = error[x]["TP"]
        tn = error[x]["TN"]
        fp = error[x]["FP"]
        fn = error[x]["FN"]
        line += str(",%d,%d,%d,%d" % (tp, tn, fp, fn))
    flip = error["AD"]["FLIP"]
    line += str(",%d" % flip)
    sys.stdout.write(line[1:])


def main(args):
    with open(args.true_mat, 'r') as fin:
        true_mat = read_matrix(fin)

    with open(args.esti_mat, 'r') as fin:
        esti_mat = read_matrix(fin)

    calculate_error(true_mat, esti_mat)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--true_mat", type=str,
                        help="True Matrix (CFMatrix format)",
                        required=True)

    parser.add_argument("-e", "--esti_mat", type=str,
                        help="Estimated Matrix (CFMatrix format)",
                        required=True)

    main(parser.parse_args())
