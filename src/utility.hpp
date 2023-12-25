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

const std::string help_info = "TO BE UPDATED!";

const std::string help_info_bak = 
"=================================== TREE-QMC ===================================\n"
"This is version 2.0.0 of TREe Embedded Quartet Max Cut (TREE-QMC).\n\n" 
"USAGE:\n"
"./TREE-QMC (-i|--input) <input file> [(-o|--output) <output file>]\n"
"           [(--polyseed) <integer>] [(--maxcutseed) <integer>]\n"
"           [(-n|--normalize) <normalization scheme>]\n"
"           [(-x|--execution) <execution mode>]\n"
"           [(-v|--verbose) <verbose mode>] [-h|--help]\n\n"
"OPTIONS:\n"
"[-h|--help]\n"
"        Prints this help message.\n"
"(-i|--input) <input file>\n"
"        Name of file containing gene trees in newick format (required)\n"
"        IMPORTANT: current implementation of TREE-QMC requires that the input\n"
"        gene trees are unrooted and binary. Thus, TREE-QMC suppresses roots\n" 
"        and randomly refines polytomies during a preprocessing phase; the\n"
"        resulting trees are written to \"<input file>.refined\".\n"
"[(-o|--output) <output file>]\n"
"        Name of file for writing output species tree (default: stdout)\n"
"[(--polyseed) <integer>]\n"
"        Seeds random number generator with <integer> prior to arbitrarily\n"
"        resolving polytomies. If <integer> is set to -1, system time is used;\n" 
"        otherwise, <integer> should be positive (default: 12345).\n"
"[(--maxcutseed) <integer>]\n"
"        Seeds random number generator with <integer> prior to calling the max\n"
"        cut heuristic but after the preprocessing phase. If <integer> is set to\n"
"        -1, system time is used; otherwise, <integer> should be positive\n"
"        (default: 1).\n"
"[(-n|--normalize) <normalization scheme>]\n"
"        Initially, each quartet is weighted by the number of input gene\n"
"        trees that induce it. At each step in the divide phase of wQMC and\n"
"        TREE-QMC, the input quartets are modified with artificial taxa. We\n"
"        introduce two normalization schemes for artificial taxa and find\n"
"        that they improve empirical performance of TREE-QMC in a simulation\n"
"        study. The best scheme is run by default. See paper for details.\n"
"        -n 0: none\n"
"        -n 1: uniform\n"
"        -n 2: non-uniform (default)\n"
"[(-x|--execution) <execution mode>]\n"
"        TREE-QMC uses an efficient algorithm that operates directly on the\n"
"        input gene trees by default. The naive algorithm, which operates on a\n"
"        set of quartets weighted based on the input gene trees, is also\n"
"        implemented for testing purposes.\n"
"        -x 0: run efficient algorithm (default)\n"
"        -x 1: run naive algorithm\n"
"[(-v|--verbose) <verbose mode>]\n"
"        -v 0: write no subproblem information (default)\n"
"        -v 1: write CSV with subproblem information (subproblem ID, parent\n"
"              problem ID, depth of recursion, total # of taxa, # of artifical\n"
"              taxa, species names)\n"
"        -v 2: write CSV with subproblem information (info from v1 plus # of\n" 
"              of elements in f, # of pruned elements in f, # of zeroes in f)\n"
"[--shared <use shared taxon data structure to normalize quartet weights>]\n"
"        Do NOT use unless there are no missing data!!!\n\n"
"Contact: Post issue to Github (https://github.com/molloy-lab/TREE-QMC/)\n"
"        or email Yunheng Han (yhhan@umd.edu) & Erin Molloy (ekmolloy@umd.edu)\n\n"
"If you use TREE-QMC in your work, please cite:\n"
"  Han and Molloy, 2022, \"TREE-QMC: Improving quartet graph construction for\n"
"  scalable and accurate species tree estimation from gene trees,\" bioRxiv,\n"
"  https://doi.org/10.1101/2022.06.25.497608.\n"
"================================================================================\n\n";

#endif
