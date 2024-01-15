TREE-QMC
========

TREE-QMC is a quartet-based method for estimating species trees from gene trees, like the popular methods [ASTRAL](https://doi.org/10.1186/s12859-018-2129-y) and [Weighted ASTRAL/ASTER](https://doi.org/10.1093/molbev/msac215). To learn more about TREE-QMC, check out [Han & Molloy, *Genome Res*, 2023](http:doi.org/10.1101/gr.277629.122).

Acknowledgements
----------------
TREE-QMC is based on the Quartet Max Cut (QMC) framework introduced by Sagi Snir and Satish Rao; see [Snir & Rao, *IEEE/ACM TCBB*, 2010](http:doi.org/10.1109/TCBB.2008.133) and [Avni, Cohen & Snir, *Syst Biol*, 2015](http:doi.org/10.1093/sysbio/syu087).

TREE-QMC implements efficient (and brute force) algorithms for the **quartet weighting schemes** introduced by Chao Zhang and Siavash Mirarab; see [Zhang & Mirarab, *Mol Biol Evol*, 2022](https://doi.org/10.1093/molbev/msac215).

[MQLib](https://github.com/MQLib/MQLib) is utilized for its max cut heuristic; see [Dunning, Gupta, & Silberholz, *INFORMS Journal on Computing*, 2018](https://doi.org/10.1287/ijoc.2017.0798).

[toms743](https://people.sc.fsu.edu/~jburkardt/cpp_src/toms743/toms743.html) is utilized for its Lambert's W approximation; see [Fritsch, Shafer & Crowley, *Communications of the ACM*, 1973](https://doi.org/10.1145/361952.361970) and [Barry, Barry & Culligan-Hensley, *ACM Transactions on Mathematical Software*, 1995](https://doi.org/10.1145/203082.203088).


Tutorial
--------
We recommend working through [this tutorial](tutorial/README.md).

Quick Tips
----------
If you are having strange problems, try removing the `\r` characters from the file before using TREE-QMC:
```
cat file.txt | tr -d '\r' > newfile.txt
```

Add the path to treeqmc to your bash profile for convenience.
```
export PATH=$PATH:"<path to treeqmc>"
```


Build
-----
To build wTREE-QMC, use commands:
```
git clone https://github.com/molloy-lab/weighted-TREE-QMC.git
cd weighted-TREE-QMC

# HANDLE MQLIB
cd external/MQLib
make
cd ../..

# TO TEST LAMBERT's W
#cd external/toms743
#g++ -std=c++11 -O2 -o toms743_test *cpp -lm
#./toms743_test
#cd ../..

g++ -std=c++11 -O2 \
    -I external/MQLib/include -I external/toms743 \
    -o treeqmc \
    src/*.cpp external/toms743/toms743.cpp \
    external/MQLib/bin/MQLib.a -lm 
```

Usage
-----
To run TREE-QMC, use command:
```
./treeqmc -i <input file> -o <output file name>
```

To see the TREE-QMC usage options, use command:
```
./treeqmc -h
```

The output help message should be
```
=================================== wTREE-QMC ===================================
USAGE:
./treeqmc (-i|--input) <input file> [(-o|--output) <output file>]
            [(-w|-weight) <weighting scheme>]
            [(-r|--support_range) <min> <max>]
            [(-c|--contract) <threshold>]
            [(-n|--normalize) <normalization scheme>]
            [(-x|--execution) <execution mode>]
            [(-v|--verbose) <verbose mode>]
            [-h|--help]

OPTIONS:
[-h|--help]
        Prints this help message.
(-i|--input) <input file>
        Name of file containing gene trees in newick format (required)
[(-o|--output) <output file>]
        Name of file for writing output species tree (default: stdout)
[(-w|--weight) <weighting scheme>]
        Weighting scheme for quartets; see paper for details
        -w n: none (default)
        -w s: support only
        -w h: hybrid of support and length
        -w l: length only
        -w f: none fast
              Refines polytomies arbitrarily so faster algorithm can be used
[(-r|--support_range) <min> <max>]
        Specifies minimum and maximum branch support values (*required* when
        using -w s or -w h options)
        + For probability or likelihood support, use: -s 0 1
        + For bootstrap support, use: -s 0 100
        + For local Bayesian support, use: -s 0.333 1 (abayes is recommended)
[(-c|--contract) <threshold>]
        Run unweighted method (-w n) after contracting internal branches with
        support less than <threshold>
[(-n|--normalize) <normalization scheme>]
        Normalization scheme for artificial taxa; see paper for details
        -n 0: none
        -n 1: uniform
        -n 2: non-uniform (default)
[(-x|--execution) <execution mode>]
        -x 0: run efficient algorithm (default)
        -x 1: run brute force algorithm for testing
        -x 2: compute weighted quartets, then exit
        -x 3: compute good and bad edges, then exit
[(-v|--verbose) <verbose mode>]
        -v 0: write no subproblem information (default)
        -v 1: write CSV with subproblem information (subproblem ID, parent
              problem ID, depth of recursion, total # of taxa, # of artifical
              taxa, species names)
        -v 2: write CSV with subproblem information (info from v1 plus # of
              of elements in f, # of pruned elements in f, # of zeroes in f)

OTHER OPTIONS:
[(--polyseed) <integer>]
        Seeds random number generator with prior to arbitrarily resolving
        polytomies. If seed is set to -1, system time is used;
        otherwise, seed should be positive (default: 12345).
[(--maxcutseed) <integer>]
        Seeds random number generator prior to calling the max cut heuristic
        but after the preprocessing phase. If seed is set to -1, system time
        is used; otherwise, seed should be positive (default: 1).
[--shared <use shared taxon data structure to normalize quartet weights>]
        Do NOT use unless there are NO missing data!!!

Contact: Post issue to Github (https://github.com/molloy-lab/weighted-TREE-QMC/)
        or email Yunheng Han (yhhan@umd.edu) & Erin Molloy (ekmolloy@umd.edu)

If you use wTREE-QMC in your work, please cite:
  Han and Molloy, 2024, https://github.com/molloy-lab/weighted-TREE-QMC.

  and

  Han and Molloy, 2023, Improving quartet graph construction for scalable
  and accurate species tree estimation from gene trees, Genome Research,
  http:doi.org/10.1101/gr.277629.122.
=================================================================================
```
