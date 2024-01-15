TREE-QMC
========

TREE-QMC is a quartet-based method for estimating species trees directly from gene trees or characters, like the popular methods [ASTRAL](https://doi.org/10.1186/s12859-018-2129-y) and [Weighted ASTRAL/ASTER](https://doi.org/10.1093/molbev/msac215). 

TREE-QMC uses a different algorithmic approach than ASTRAL/ASTER. To learn more about the TREE-QMC algorithm, check out [Han & Molloy, *Genome Res*, 2023](http:doi.org/10.1101/gr.277629.122).

Acknowledgements
----------------
TREE-QMC is based on the Quartet Max Cut (QMC) framework introduced by Sagi Snir and Satish Rao; see [Snir & Rao, *IEEE/ACM TCBB*, 2010](http:doi.org/10.1109/TCBB.2008.133) and [Avni, Cohen & Snir, *Syst Biol*, 2015](http:doi.org/10.1093/sysbio/syu087). TREE-QMC contributes fast algorithms for constructing the quartet graph directly from the input trees, rather than explicitly enumerating all quartets or sampling quartets.

TREE-QMC now implements efficient (and brute force) algorithms for the **quartet weighting schemes** introduced by Chao Zhang and Siavash Mirarab; see [Zhang & Mirarab, *Mol Biol Evol*, 2022](https://doi.org/10.1093/molbev/msac215).

TREE-QMC uses [MQLib](https://github.com/MQLib/MQLib) for its max cut heuristic; see [Dunning, Gupta, & Silberholz, *INFORMS Journal on Computing*, 2018](https://doi.org/10.1287/ijoc.2017.0798).

TREE-QMC uses [toms743](https://people.sc.fsu.edu/~jburkardt/cpp_src/toms743/toms743.html) for its Lambert's W approximation; see [Fritsch, Shafer & Crowley, *Communications of the ACM*, 1973](https://doi.org/10.1145/361952.361970) and [Barry, Barry & Culligan-Hensley, *ACM Transactions on Mathematical Software*, 1995](https://doi.org/10.1145/203082.203088).

Build
-----
To build TREE-QMC, use commands:
```
git clone https://github.com/molloy-lab/weighted-TREE-QMC.git
cd weighted-TREE-QMC/external/MQLib
make
cd ../..
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
./treeqmc -i <input file>
```

Tips
----
If you are having strange problems, try removing the `\r` characters from the input files and trying again:
```
cat file.txt | tr -d '\r' > newfile.txt
```
For convenience, add treeqmc to your shell profile. For bash, open `~/.bash_profile` and add
```
export PATH=$PATH:"<path to treeqmc>"
```

Tutorial
--------
Check out [this tutorial](tutorial/README.md).

Options
-------

To see the TREE-QMC usage options, use command:
```
./treeqmc -h
```

The output help message should be
```
=================================== TREE-QMC ===================================
BASIC USAGE:
./treeqmc (-i|--input) <input file>

Help Options:
[-h|--help]
        Prints this help message.

Input Options:
(-i|--input) <input file>
        File with gene trees in newick format (required)
[(--chars)]
        Input are characters rather than trees
        Missing states are N, -, and ?
[(--bp)]
        Input are binary characters i.e. bipartitions
[(-a|-mapping) <mapping file>]
        File with individual/leaf names (1st col) mapped to species (2nd col)
[(--root) <list of species separated by commas>]
        Root species tree at given species if possible
[(--rootonly) <species tree file>]
        Root species tree in file and then exit

[(--supportonly) <species tree file>]
        Compute quartet support for species tree in file and then exit

Output Options:
[(-o|--output) <output file>]
        File for writing output species tree (default: stdout)
[(--support)]
        Compute quartet support for output species tree
[(--writetable) <table file>]
        Write branch and quartet support information to CSV
[(--char2tree)]
        Write character matrix as trees (newick strings) to output and exit

Algorithm Options:
[(--hybrid)]
        Use hybrid weighting scheme (-w h)
[(--fast)]
        Use fast algorithm that does not support weights or polytomies (-w f)
[(-B|--bayes)]
       Use presets for bayesian support (-n 0.333 -x 1.0 -d 0.333)
[(-L|--lrt)]
       Use presets for likelihood support (-n 0.0 -x 1.0 -d 0.0)
[(-S|--bootstrap)]
       Use presets for boostrap support (-n 0 -x 100 -d 0)
[(-c|--contract) <float>]
       Contract internal branches with support less than specified threshold

       after mapping suport to the interval 0 to 1Advanced Options:
[(-w|--weight) <character>]
        Weighting scheme for quartets; see paper for details
        -w n: none (default)
        -w h: hybrid of support and length (recommended)
        -w s: support only
        -w l: length only
        -w f: none fast
              Refines polytomies arbitrarily so faster algorithm can be used
[(-n|--min) <float>]
        Minimum value of input branch support (default: 0.0)
[(-x|--max) <float>]
        Maximum value of input branch support (default: 1.0)
[(-d|--default) <float>]
        Default value of input branch support when not provided (default: 0.0)
[(--norm_atax) <integer>]
        Normalization scheme for artificial taxa; see paper for details
        -n 0: none
        -n 1: uniform
        -n 2: non-uniform (default)
[(-e|--execution) <execution mode>]
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
================================================================================
```
