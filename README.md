<a href="url"><img src="https://github.com/molloy-lab/TREE-QMC/blob/main/external/logo.png" align="left" height="135" width="135" ></a>

TREE-QMC
========

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/tree-qmc/README.html)

TREE-QMC is a quartet-based method for estimating species trees directly from gene trees or characters, like the popular method ASTRAL methods; see [Han & Molloy, 2023](https://doi.org/10.1101/gr.277629.122), [2025](https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syaf009/8042591?utm_source=authortollfreelink&utm_campaign=sysbio&utm_medium=email&guestAccessKey=24c3e656-5c43-482a-8cdd-c7fa757798d0). Unlike ASTRAL, TREE-QMC uses a different algorithmic approach, based on the Quartet Max Cut (QMC) framework of Snir and Rao. This approach is particularly beneficial for phylogenomic data sets with high **missingness**. Additionally, TREE-QMC can be used to reconstruct a **tree of blobs** and evaluate signals of gene flow and other network-level evolutionary processes; see [Dai et al., 2024](https://doi.org/10.1101/2025.11.05.686850).


TUTORIALS
---------
Check out:
+ [tutorial for tree of blobs (TOB) reconstruction and exploration of network-like evolution](tutorial/tree-of-blobs/README.md)
+ [tutorial for species/population tree estimation from gene trees](tutorial/gene-trees/README.md)
+ [tutorial for species/population tree estimation from multi-labeled gene trees](tutorial/multi-gene-trees/README.md)
+ [tutorial for species/population tree estimation character matrices](tutorial/characters/README.md)
+ [tutorial for Partitioned Coalescence Support (PCS)](tutorial/characters/README.md)
+ [tutorial for quartet inputs](tutorial/quartets/README.md)


BUILD
-----
Check out the [tree of blobs tutorial](tutorial/tree-of-blobs/README.md) for build instructions.
To build TREE-QMC *without* tree of blobs functionality use the following commands:
```
git clone https://github.com/molloy-lab/TREE-QMC
cd TREE-QMC/external/MQLib
make
cd ../../ && mkdir -p build && cd build
g++ -std=c++11 -O2 \
    -I ../external/MQLib/include \
    -I ../external/toms743 \
    -o tree-qmc \
    ../src/*.cpp \
    ../external/toms743/toms743.cpp \
    ../external/MQLib/bin/MQLib.a \
    -lm \
    -DVERSION=\"$(cat ../version.txt)\"
```
Lastly, add `build` directory to your `$PATH` environment variable by typing
```
export PATH="$(pwd):$PATH"
```
so that your system can find `tree-qmc` from other directories.
Even better, add the following line
```
export PATH="<path to TREE-QMC>:$PATH"
```
to your `~/.bash_profile` file or your `~/.bashrc` file so that your system can find `tree-qmc` whenever you start a new terminal instance.


USAGE
-----
To see the TREE-QMC usage options, type command:
```
tree-qmc -h
```

The output help message should be
```
=================================== TREE-QMC ===================================
USAGE:
tree-qmc -i <input gene trees> -o <output species tree>

Note: If directory containing tree-qmc is NOT on your $PATH, replace
      tree-qmc with <path to tree-qmc>/tree-qmc in command above

Help Options:
[-h|--help]
        Prints this help message.

Input Options:
(-i|--input) <input file>
        Input file, typically with gene trees in newick format (required)
[(--quartets)]
        Input are (weighted) quartets, either in qCF format or format below
        Note: only n0 and n2 normalization are implemented for quartet input
[(--quartetformat) <format string>]
        Specify a format of the input quartets
        Examples: "((___,___),(___,___));___" for ((A,B),(C,D));1.234 (default)
                  or "___,___|___,___:___" for A,B|C,D:1.234
[(--chars)]
        Input are characters in fasta, phylip, or nexus format
        Missing states are N, -, and ?
[(--bp)]
        Input are binary (0/1) characters in fasta, phylip, or nexus format
        Missing states are N, -, and ?
[(-a|-mapping) <mapping file>]
        File with individual/leaf names mapped to species names
        Example: ind1 taxA
                 ind2 taxA
                 ind3 taxB
                 ...

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
        Compute the tree of blobs from the input gene trees
[(--store_pvalue)]
        Store min p-value found for each branch, then exit
[(--3f1a)]
        Use 3-fix-1-alter algorithm for minimum p-value search
[(--iter_limit_blob) <non-negative integer>]
        Maximum number of iterations for default (bipartition) search algorithm for
        min p-value (default: two times the number of taxa squared)
        Set to 0 to perform exhaustive search for min p-value
[(--load_pvalue)]
        Load tree with p-values and contract branches based on alpha, beta settings
[(--alpha <float number>)]
        Hyperparameter for hypothesis testing with tree-test (default: 1e-7)
[(--beta <float number>)]
        Hyperparameter for hypothesis testing with star-test (default: 0.95)
[(--blobsearchonly) <base tree file>]
        Perform hypothesis testing on input base tree

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


TIPS
----
Hidden symbols in your input data file can cause very strange problems!! If you are having strange problems, try removing the hidden `\r` symbol from your input files with the following command 
```
cat <input file> | tr -d '\r' > <clean input file>
```
and then running TREE-QMC on the cleaned data file. This is a common issue if your input files were created on a Windows operating system.


ACKNOWLEDGEMENTS
----------------
+ TREE-QMC is based on the Quartet Max Cut (QMC) framework introduced by Sagi Snir and Satish Rao; see [Snir & Rao, *IEEE/ACM TCBB*, 2010](http:doi.org/10.1109/TCBB.2008.133) and [Avni, Cohen & Snir, *Syst Biol*, 2015](http:doi.org/10.1093/sysbio/syu087). It contributes fast algorithms for constructing the quartet graph directly from the input trees, rather than enumerating quartets or sampling quartets.

+ TREE-QMC enables quartet branch support calculations with the approach of [Sayyari & Mirarab, *Mol Biol Evol*, 2016](https://doi.org/10.1093/molbev/msw079).

+ Weighted TREE-QMC leverages the **quartet weighting schemes** introduced by Chao Zhang and Siavash Mirarab; see [Zhang & Mirarab, *Mol Biol Evol*, 2022](https://doi.org/10.1093/molbev/msac215). 

+ TOB-QMC uses the hypothesis tests implemented in [TINNiK](https://cran.r-project.org/web/packages/MSCquartets/vignettes/TINNIK.html) to infer tree of blobs; see [Allman, et al., *Algorithms Mol Biol*, 2024](https://doi.org/10.1186/s13015-024-00266-2)

+ TREE-QMC uses [MQLib](https://github.com/MQLib/MQLib) for its max cut heuristic; see [Dunning, Gupta, & Silberholz, *INFORMS J Computing*, 2018](https://doi.org/10.1287/ijoc.2017.0798).

+ TREE-QMC uses [toms743](https://people.sc.fsu.edu/~jburkardt/cpp_src/toms743/toms743.html) for its Lambert's W approximation; see [Fritsch, Shafer & Crowley, *Communications of the ACM*, 1973](https://doi.org/10.1145/361952.361970) and [Barry, Barry & Culligan-Hensley, *ACM Trans Math Software*, 1995](https://doi.org/10.1145/203082.203088).
