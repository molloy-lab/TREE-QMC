#ifndef INSTANCE_HPP
#define INSTANCE_HPP

#include "utility.hpp"
#include "tree.hpp"
#include "charmat.hpp"

class Instance {
    public:
        Instance(int argc, char **argv);
        ~Instance();
        long long solve();
        SpeciesTree *get_solution();
        void output_solution();
        std::string get_execution_mode();
    private:
        std::vector<Tree *> input;
        Dict *dict;
        SpeciesTree *output;
        std::unordered_set<std::string> outgroup_taxon_set;
        std::unordered_map<std::string, std::string> indiv2taxon;
        std::string input_file, output_file, mapping_file, stree_file, table_file;
        std::string root_str;
        std::string normal_mode, weight_mode, execute_mode, taxa_mode, score_mode, data_mode, brln_mode;
        unsigned long int refine_seed, cut_seed, iter_limit;
        weight_t support_low, support_high, support_default, support_threshold;
        bool contract, char2tree, rootonly;
        int parse(int argc, char **argv);
        void prepare_root_taxa();
        void prepare_indiv2taxon_map();
        void input_trees();
        void input_matrix();
        void prepare_trees();
        void refine_trees();
};

extern std::ofstream subproblem_csv;
extern std::string verbose;
extern unsigned long long count[8];

#endif
