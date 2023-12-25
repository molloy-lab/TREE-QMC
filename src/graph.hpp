#ifndef GRAPH_HPP
#define GRAPH_HPP

#include "utility.hpp"
#include "dict.hpp"
#include "tree.hpp"
#include "heuristics/maxcut/burer2002.h"
#include "problem/instance.h"
#include "problem/max_cut_instance.h"

class Graph {
    public:
        Graph(std::vector<Tree *> trees, Taxa &subset);
        Graph(std::vector<Tree *> trees, Taxa &subset, std::string weighting);
        Graph(std::unordered_map<quartet_t, weight_t> &quartets, Taxa &subset);
        ~Graph();
        std::string to_string();
        weight_t get_cut(std::vector<index_t> *A, std::vector<index_t> *B, unsigned long int iter_limit);
    private:
        index_t size;
        std::unordered_map<index_t, index_t> index2index;
        std::vector<index_t> indices;
        weight_t ***graph;
        weight_t sdp_cut(weight_t alpha, std::vector<index_t> *A, std::vector<index_t> *B, unsigned long int iter_limit);
};

class QuartetGraphMaxCutCallback : public MaxCutCallback {
    public:
        QuartetGraphMaxCutCallback(unsigned long int iter_limit);
        ~QuartetGraphMaxCutCallback();
        bool Report(const MaxCutSimpleSolution& solution, bool newBest, double runtime);
        bool Report(const MaxCutSimpleSolution& solution, bool newBest, double runtime, int iter);
    private:
        unsigned long int iter_limit;
};

extern std::ofstream subproblem_csv;
extern std::string verbose;
extern unsigned long long count[8];

#endif
