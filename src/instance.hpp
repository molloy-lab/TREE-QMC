#ifndef INSTANCE_HPP
#define INSTANCE_HPP

#include "utility.hpp"
#include "tree.hpp"

class Instance {
    public:
        Instance(int argc, char **argv);
        ~Instance();
        long long solve();
        SpeciesTree *get_solution();
        void output_solution();
        std::string get_execution_mode();
    private:
        std::string input_file, output_file, normal, weight, execute, taxa_mode, output_quartets, score_mode;
        unsigned long int refine_seed, cut_seed, trc, iter_limit;
        weight_t support_low, support_high, threshold;
        bool contract;
        std::vector<Tree *> input;
        Dict *dict;
        SpeciesTree *output;
        bool parse(int argc, char **argv);
        std::size_t input_trees();
        void resolve_trees();
        void prepare_trees();
};

extern std::ofstream subproblem_csv, quartets_txt, good_edges_txt, bad_edges_txt;
extern std::string verbose;
extern unsigned long long count[8];

#endif
