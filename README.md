<a href="url"><img src="https://github.com/molloy-lab/TREE-QMC/blob/main/external/logo.png" align="left" height="135" width="135" ></a>

TREE-QMC
========

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/tree-qmc/README.html)

TREE-QMC is a quartet-based method for estimating species trees directly from gene trees or characters, like the popular methods [ASTRAL](https://doi.org/10.1186/s12859-018-2129-y) and [Weighted ASTRAL/ASTER](https://doi.org/10.1093/molbev/msac215). However, TREE-QMC uses a different algorithmic approach than ASTRAL/ASTER, based on the Quartet Max Cut (QMC) framework of Snir and Rao, that is particularly beneficial for large phylogenomic data sets with high missingness.


Tutorials
---------
Check out: 
+ [tutorial for gene trees](tutorial/gene-trees/README.md)
+ [tutorial for multi-labeled gene trees](tutorial/multi-gene-trees/README.md)
+ [tutorial for character matrices and Partitioned Coalescence Support (PCS) analyses](tutorial/characters/README.md)
+ [tutorial for tree of blobs (TOB)](tutorial/tree-of-blobs/README.md)

Tips
----
**Tip #1.** Add the directory containing the `tree-qmc` binary to your `PATH` variable so that you can type `tree-qmc` instead of `<path to tree-qmc binary>/tree-qmc`. For bash, this can be done by adding
```
export PATH=$PATH:"path to tree-qmc binary"
```
to the `~/.bash_profile` file or the `~/.bashrc` file.

**Tip #2.** Hidden symbols in your input data file can cause strange problems. Try removing `\r` symbols from the input files with the following command 
```
cat <input data> | tr -d '\r' > <clean input data>
```
and then trying to run TREE-QMC again.

BUILD
-----
We recommend you install TREE-QMC with bioconda or download an executable for convenience. To build TREE-QMC (without tree of blobs functionality) use the following commands:
```
git clone htts://github.com/molloy-lab/TREE-QMC
cd TREE-QMC
cd external/MQLib
make
cd ../../
g++ -std=c++11 -O2 \
    -I external/MQLib/include \
    -I external/toms743 \
    -o tree-qmc \
    src/*.cpp \
    external/toms743/toms743.cpp \
    external/MQLib/bin/MQLib.a -lm \
    -DVERSION=\"$(cat version.txt)\"
```
Otherwise see the tree of blobs tutorial for build instructions.

Usage
-----
To see the TREE-QMC usage options, use command:
```
./tree-qmc -h
```

The output help message should be
```
=================================== TREE-QMC ===================================
SPECIES TREE BASIC USAGE:
tree-qmc -i <input gene trees> -o <output species tree>

TREE OF BLOBS (TOB) BASIC USAGE:
tree-qmc --store_pvalue -i <input gene trees> -o <output tree with p-values>
tree-qmc --blob --load_pvalue -i <input tree with p-values> -o <output TOB>

**If the directory containing the tree-qmc binary is not part of $PATH, replace
  tree-qmc with <path to tree-qmc binary>/tree-qmc in the command above**

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
        Missing states are N, -, and ?
[(-a|-mapping) <mapping file>]
        File with individual/leaf names (1st col) mapped to species (2nd col)

Output Options:
[(-o|--output) <output file>]
        File for writing output species tree (default: stdout)
[(--override)]
        Override output file if it already exists
[(--root) <list of species separated by commas>]
        Root species tree at given species if possible
[(--support)]
        Compute quartet support for output species tree
[(--writetable) <table file>]
        Write branch and quartet support information to CSV

Weighted Algorithm Options:
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
       after mapping suport to the interval 0 to 1

Branch Support and Utility Options:
[(--char2tree)]
        Write character matrix as trees (newick strings) to output and exit
[(--rootonly) <species tree file>]
        Root species tree in file and then exit
[(--supportonly) <species tree file>]
        Compute quartet support for species tree in file and then exit
[(--pcsonly) <species tree file>]
        Compute partitioned coalescent support (PCS) for specified branch in
        species tree in file (anotate branch with PCS) and then exit

Tree of Blobs Options:
[(--blob)]
        Compute the tree of blobs directed from the input gene trees
[(--store_pvalue)]
        Run TREE-QMC and store min p-value found for each branch, then exit
[(--3f1a)]
        Use 3-fix-1-alter search algorithm for minimum p-value
[(--iter_limit_blob) <non-negative integer>]
        Maximum number of iterations for default bipartition search algorithm for
        min p-value (default: two times the number of taxa squared)
        Set to 0 to perform exhaustive search for min p-value
[(--load_pvalue)]
        Load tree with p-values and contract branches based on alpha, beta settings
[(--alpha <float number>)]\n"
        Hyperparameter for hypothesis testing with tree-test (default: 1e-7)
[(--beta <float number>)]\n"
        Hyperparameter for hypothesis testing with star-test (default: 0.95)

Experimental/Advanced Options:
[(-w|--weight) <character>]
        Weighting scheme for quartets; see paper for details
        -w n: none (default)
        -w h: hybrid of support and length (recommended)
        -w s: support only
        -w l: length only
        -w f: none fast
              Refines polytomies arbitrarily so faster algorithm can be used
[(-n|--min) <float>]
        Minimum value of input branch support
[(-x|--max) <float>]
        Maximum value of input branch support
[(-d|--default) <float>]
        Default branch support to use if branch input tree is missing support
        Default branch length is 0.0
[(--norm_atax) <integer>]
        Normalization scheme for artificial taxa; see paper for details
        --norm_atax 0: none
        --norm_atax 1: uniform
        --norm_atax 2: non-uniform (default)
[(-e|--execution) <execution mode>]
        -e 0: run efficient algorithm (default)
        -e 1: run brute force algorithm for testing
        -e 2: compute weighted quartets, then exit
        -e 3: compute good and bad edges, then exit
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


Contact: Post issues to Github at https://github.com/molloy-lab/TREE-QMC/


If you use TREE-QMC, please cite:

  Han and Molloy, 2023, Improving quartet graph construction for scalable
  and accurate species tree estimation from gene trees, Genome Res,
  http:doi.org/10.1101/gr.277629.122

If you use weighted TREE-QMC in your work, please cite:

  Han and Molloy, 2025, Improved robustness to gene tree incompleteness,
  estimation errors, and systematic homology errors with weighted TREE-QMC,
  Syst Biol, https://doi.org/10.1093/sysbio/syaf009

If you use TOB-QMC in your work, please cite:

  Dai, Han, Molloy, 2025, Fast and accurate tree of blobs reconstruction under
  the network multispcies coalescent, bioRxiv,
  https://doi.org/10.1101/2025.11.05.686850

  Allman et al., 2024, TINNiK: inference of the tree of blobs of a species
  network under the coalescent model, Algorithms Mol Biol,
  https://doi.org/10.1186/s13015-024-00266-2
================================================================================
```

Acknowledgements
----------------
+ TREE-QMC is based on the Quartet Max Cut (QMC) framework introduced by Sagi Snir and Satish Rao; see [Snir & Rao, *IEEE/ACM TCBB*, 2010](http:doi.org/10.1109/TCBB.2008.133) and [Avni, Cohen & Snir, *Syst Biol*, 2015](http:doi.org/10.1093/sysbio/syu087). It contributes fast algorithms for constructing the quartet graph directly from the input trees, rather than enumerating quartets or sampling quartets.

+ TREE-QMC enables quartet branch support calculations with the approach of [Sayyari & Mirarab, *Mol Biol Evol*, 2016](https://doi.org/10.1093/molbev/msw079).

+ Weighted TREE-QMC leverages the **quartet weighting schemes** introduced by Chao Zhang and Siavash Mirarab; see [Zhang & Mirarab, *Mol Biol Evol*, 2022](https://doi.org/10.1093/molbev/msac215). 

+ TREE-QMC uses [MQLib](https://github.com/MQLib/MQLib) for its max cut heuristic; see [Dunning, Gupta, & Silberholz, *INFORMS J Computing*, 2018](https://doi.org/10.1287/ijoc.2017.0798).

+ TREE-QMC uses [toms743](https://people.sc.fsu.edu/~jburkardt/cpp_src/toms743/toms743.html) for its Lambert's W approximation; see [Fritsch, Shafer & Crowley, *Communications of the ACM*, 1973](https://doi.org/10.1145/361952.361970) and [Barry, Barry & Culligan-Hensley, *ACM Trans Math Software*, 1995](https://doi.org/10.1145/203082.203088).

+ TOB-QMC uses the hypothesis tests implemented in [TINNiK](https://cran.r-project.org/web/packages/MSCquartets/vignettes/TINNIK.html) to infer tree of blobs; see [Allman, et al., *Algorithms Mol Biol*, 2024](https://doi.org/10.1186/s13015-024-00266-2)

