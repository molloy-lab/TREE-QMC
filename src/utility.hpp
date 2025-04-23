#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iomanip>
#include <queue>
#include <stack>
#include <deque>
#include <RInside.h>

#define INDEX_WIDTH 65536
#define ERROR_BOUND 1e-6

typedef int16_t index_t;
typedef double weight_t;
typedef uint64_t quartet_t;

class Matrix {
    public:
        static weight_t **new_mat(index_t size);
        static void delete_mat(weight_t **m, index_t size);
        static std::string display_mat(weight_t **m, index_t size);
        static weight_t diff_mat(weight_t **m1, weight_t **m2, index_t size);
};

quartet_t join(index_t *quartet);
index_t *split(quartet_t quartet);
bool s2d(std::string s, weight_t *r);
bool s2ul(std::string s, unsigned long int *r);
weight_t *init(index_t size);
weight_t pvalue(index_t *indices);
std::vector<weight_t> pvalue_all(index_t *indices);
weight_t pvalue(weight_t *qCF);

const std::string help_info = 
"=================================== TREE-QMC ===================================\n"
"BASIC USAGE:\n"
"tree-qmc (-i|--input) <input file>\n\n"
"**If the directory containing the tree-qmc binary is not part of $PATH, replace\n"
"eplace tree-qmc with <path to tree-qmc binary>/tree-qmc in the command above**\n\n"
"Help Options:\n"
"[-h|--help]\n"
"        Prints this help message.\n\n"
"Input Options:\n"
"(-i|--input) <input file>\n"
"        File with gene trees in newick format (required)\n"
"[(--chars)]\n"
"        Input are characters rather than trees\n"
"        Missing states are N, -, and ?\n"
"[(--bp)]\n"
"        Input are binary characters i.e. bipartitions\n"
"        Missing states are N, -, and ?\n"
"[(-a|-mapping) <mapping file>]\n"
"        File with individual/leaf names (1st col) mapped to species (2nd col)\n"
"[(--root) <list of species separated by commas>]\n"
"        Root species tree at given species if possible\n"
"[(--rootonly) <species tree file>]\n"
"        Root species tree in file and then exit\n"
"[(--supportonly) <species tree file>]\n"
"        Compute quartet support for species tree in file and then exit\n\n"
"[(--pcsonly) <species tree file>]\n"
"        Compute partitioned coalescent support (PCS) for specified branch in\n"
"        species tree in file (anotate branch with PCS) and then exit\n\n"
"Output Options:\n"
"[(-o|--output) <output file>]\n"
"        File for writing output species tree (default: stdout)\n"
"[(--support)]\n"
"        Compute quartet support for output species tree\n"
"[(--writetable) <table file>]\n"
"        Write branch and quartet support information to CSV\n"
"[(--char2tree)]\n"
"        Write character matrix as trees (newick strings) to output and exit\n\n"
"Algorithm Options:\n"
"[(--hybrid)]\n"
"        Use hybrid weighting scheme (-w h)\n"
"[(--fast)]\n"
"        Use fast algorithm that does not support weights or polytomies (-w f)\n"
"[(-B|--bayes)]\n"
"       Use presets for bayesian support (-n 0.333 -x 1.0 -d 0.333)\n"
"[(-L|--lrt)]\n"
"       Use presets for likelihood support (-n 0.0 -x 1.0 -d 0.0)\n"
"[(-S|--bootstrap)]\n"
"       Use presets for boostrap support (-n 0 -x 100 -d 0)\n"
"[(-c|--contract) <float>]\n"
"       Contract internal branches with support less than specified threshold\n"
"       after mapping suport to the interval 0 to 1\n\n"
"Advanced Options:\n"
"[(-w|--weight) <character>]\n"
"        Weighting scheme for quartets; see paper for details\n"
"        -w n: none (default)\n"
"        -w h: hybrid of support and length (recommended)\n"        
"        -w s: support only\n"
"        -w l: length only\n"
"        -w f: none fast\n"
"              Refines polytomies arbitrarily so faster algorithm can be used\n"
"[(-n|--min) <float>]\n"
"        Minimum value of input branch support\n"
"[(-x|--max) <float>]\n"
"        Maximum value of input branch support\n"
"[(-d|--default) <float>]\n"
"        Default branch support to use if branch input tree is missing support\n"
"        Default branch length is 0.0\n"
"[(--norm_atax) <integer>]\n"
"        Normalization scheme for artificial taxa; see paper for details\n"
"        --norm_atax 0: none\n"
"        --norm_atax 1: uniform\n"
"        --norm_atax 2: non-uniform (default)\n"
"[(-e|--execution) <execution mode>]\n"
"        -e 0: run efficient algorithm (default)\n"
"        -e 1: run brute force algorithm for testing\n"
"        -e 2: compute weighted quartets, then exit\n"
"        -e 3: compute good and bad edges, then exit\n"
"[(-v|--verbose) <verbose mode>]\n"
"        -v 0: write no subproblem information (default)\n"
"        -v 1: write CSV with subproblem information (subproblem ID, parent\n"
"              problem ID, depth of recursion, total # of taxa, # of artifical\n"
"              taxa, species names)\n"
"        -v 2: write CSV with subproblem information (info from v1 plus # of\n" 
"              of elements in f, # of pruned elements in f, # of zeroes in f)\n"
"[(--polyseed) <integer>]\n"
"        Seeds random number generator with prior to arbitrarily resolving\n"
"        polytomies. If seed is set to -1, system time is used;\n" 
"        otherwise, seed should be positive (default: 12345).\n"
"[(--maxcutseed) <integer>]\n"
"        Seeds random number generator prior to calling the max cut heuristic\n"
"        but after the preprocessing phase. If seed is set to -1, system time\n"
"        is used; otherwise, seed should be positive (default: 1).\n"
"[--shared <use shared taxon data structure to normalize quartet weights>]\n"
"        Do NOT use unless there are NO missing data!!!\n\n"
"Contact: Post issue to Github (https://github.com/molloy-lab/TREE-QMC/)\n"
"        or email Yunheng Han (yhhan@umd.edu) & Erin Molloy (ekmolloy@umd.edu)\n\n"
"If you use TREE-QMC or weighted TREE-QMC in your work, please cite:\n\n"
"  Han and Molloy, 2023, Improving quartet graph construction for scalable\n"
"  and accurate species tree estimation from gene trees, Genome Research,\n"
"  http:doi.org/10.1101/gr.277629.122.\n\n"
"  Han and Molloy, 2024, Improved robustness to gene tree incompleteness,\n" 
"  estimation errors, and systematic homology errors with weighted TREE-QMC,\n" 
"  bioRxiv, https://doi.org/10.1101/2024.09.27.615467.\n"
"================================================================================\n\n";
#endif
