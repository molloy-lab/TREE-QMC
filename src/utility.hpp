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

const std::string help_info = 
"=================================== wTREE-QMC ===================================\n"
"USAGE:\n"
"./wTREE-QMC (-i|--input) <input file> [(-o|--output) <output file>]\n"
"            [(-w|-weight) <weighting scheme>]\n"
"            [(-r|--support_range) <min> <max>]\n"
"            [(-c|--contract) <threshold>]\n"
"            [(-n|--normalize) <normalization scheme>]\n"
"            [(-x|--execution) <execution mode>]\n"
"            [(-v|--verbose) <verbose mode>]\n"
"            [-h|--help]\n\n"
"OPTIONS:\n"
"[-h|--help]\n"
"        Prints this help message.\n"
"(-i|--input) <input file>\n"
"        Name of file containing gene trees in newick format (required)\n"
"[(-o|--output) <output file>]\n"
"        Name of file for writing output species tree (default: stdout)\n"
"[(-w|--weight) <weighting scheme>]\n"
"        Weighting scheme for quartets; see paper for details\n"
"        -w n: none (default)\n"
"        -w s: support only\n"
"        -w h: hybrid of support and length\n"
"        -w l: length only\n"
"        -w f: none fast\n"
"              Refines polytomies arbitrarily so faster algorithm can be used\n"
"[(-r|--support_range) <min> <max>]\n"
"        Specifies minimum and maximum branch support values (*required* when\n"
"        using -w s or -w h options)\n"
"        + For probability or likelihood support, use: -s 0 1\n"
"        + For bootstrap support, use: -s 0 100\n"
"        + For local Bayesian support, use: -s 0.333 1 (abayes is recommended)\n"
"[(-c|--contract) <threshold>]\n"
"        Run unweighted method (-w n) after contracting internal branches with\n"
"        support less than <threshold>\n"
"[(-n|--normalize) <normalization scheme>]\n"
"        Normalization scheme for artificial taxa; see paper for details\n"
"        -n 0: none\n"
"        -n 1: uniform\n"
"        -n 2: non-uniform (default)\n"
"[(-x|--execution) <execution mode>]\n"
"        -x 0: run efficient algorithm (default)\n"
"        -x 1: run brute force algorithm for testing\n"
"        -x 2: compute weighted quartets, then exit\n"
"        -x 3: compute good and bad edges, then exit\n"
"[(-v|--verbose) <verbose mode>]\n"
"        -v 0: write no subproblem information (default)\n"
"        -v 1: write CSV with subproblem information (subproblem ID, parent\n"
"              problem ID, depth of recursion, total # of taxa, # of artifical\n"
"              taxa, species names)\n"
"        -v 2: write CSV with subproblem information (info from v1 plus # of\n" 
"              of elements in f, # of pruned elements in f, # of zeroes in f)\n\n"
"OTHER OPTIONS:\n"
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
"Contact: Post issue to Github (https://github.com/molloy-lab/weighted-TREE-QMC/)\n"
"        or email Yunheng Han (yhhan@umd.edu) & Erin Molloy (ekmolloy@umd.edu)\n\n"
"If you use wTREE-QMC in your work, please cite:\n"
"  Han and Molloy, 2024, https://github.com/molloy-lab/weighted-TREE-QMC.\n\n"
"  and\n\n"
"  Han and Molloy, 2023, Improving quartet graph construction for scalable\n"
"  and accurate species tree estimation from gene trees, Genome Research,\n"
"  http:doi.org/10.1101/gr.277629.122.\n"
"=================================================================================\n\n";

#endif
