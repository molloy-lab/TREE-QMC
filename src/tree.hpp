#ifndef TREE_HPP
#define TREE_HPP

#include "utility.hpp"
#include "dict.hpp"
#include "taxa.hpp"

class Node {
    friend class Tree;
    friend class SpeciesTree;
    public:
        Node(index_t index);
        Node(index_t index, bool isfake);
        ~Node();
        void new_states(index_t size);
        void delete_states();
    private:
        Node *parent;
        std::vector<Node *> children;
        index_t index, size, depth;
        weight_t /* **doublet, */ *singlet, s1, s2, support, length;
        // std::map<index_t, weight_t> doublet;
        std::vector<std::pair<index_t, weight_t>> *doublet;
        static weight_t get_doublet(weight_t *singlet, weight_t s1, weight_t s2, index_t x, index_t y);
        weight_t get_doublet(index_t a, index_t b);
        void add_doublet(index_t a, index_t b, weight_t c);
        // for weighted quartets:
        weight_t length_, support_[2], plength, tdoublet[2], tdoublet_[2];
        weight_t *ssinglet, *ssinglet_, *pdoublet[2], *ptriplet[2], *mdoublet[2], *mdoublet_[2];
        weight_t **sdoublet[2], **sdoublet_[2], **striplet[2];
        bool isfake;
        weight_t f[3];
};

class Tree {
    friend class SpeciesTree;
    public:
        Tree();
        Tree(const std::string &newick, Dict *dict);
        Tree(const std::string &newick, Dict *dict, std::string support);
        virtual ~Tree();
        std::string to_string();
        std::string to_string_basic();
        size_t resolve();
        void prepare(std::string weight, weight_t low, weight_t high, bool contract, weight_t threshold);
        index_t size();
        std::unordered_map<index_t, index_t> &get_indices();
        weight_t ***build_graph(Taxa &subset);
        weight_t ***build_wgraph(Taxa &subset);
        void get_quartets(std::unordered_map<quartet_t, weight_t> *quartets);
        void get_wquartets(std::unordered_map<quartet_t, weight_t> *quartets);
        void get_wquartets_(std::unordered_map<quartet_t, weight_t> *quartets);
        std::string to_string(std::unordered_map<quartet_t, weight_t> &quartets);
        void test(Taxa &subset);
    protected:
        Node *root;
        std::unordered_map<index_t, Node*> index2node;
        Dict *dict;
        index_t pseudonym();
        std::string display_tree(Node *root);
        std::string display_tree_basic(Node *root);
        std::string display_tree_index(Node *root);
    private:
        index_t pseudonyms;
        std::unordered_map<index_t, index_t> indices;
        void clear_states(Node *root);
        void build_states(Node *root, Taxa &subset);
        void depth(Node *root, index_t depth);
        weight_t get_doublet(Node *subtree, index_t x, index_t y, bool complement);
        void sa_doublet(Node *root, weight_t sum, index_t x);
        weight_t aa_doublet(Node *root, index_t x, index_t y);
        static bool cmp(const std::pair<index_t, weight_t> &a, const std::pair<index_t, weight_t> &b);
        void sort_doublet(Node *root);
        std::unordered_set<index_t> bad_edges(Node *root, Taxa &subset, weight_t ***graph);
        void good_edges(Node *root, Taxa &subset, weight_t ***graph);
        Node *build_tree(const std::string &newick);
        Node *build_tree(const std::string &newick, std::string support);
        Node *build_subtree_from(Node *root);
        size_t resolve_tree(Node *root);
        void prepare_tree(Node *root, std::string weight, weight_t low, weight_t high, bool contract, weight_t threshold);
        void resolve_support(Node *root);
        void add_indices(Node *root, std::vector<index_t> &indices);
        void get_leaves(Node *root, std::vector<Node *> *leaves);
        void get_depth(Node *root, index_t depth);
        //for weighted quartets:
        void clear_wstates(Node *root);
        void build_wstates(Node *root, Taxa &subset);
        void build_ssinglet(Node *root, Taxa &subset);
        void build_ssinglet_(Node *root);
        void build_sdoublet(Node *root);
        void build_sdoublet_(Node *root);
        void build_striplet(Node *root);
        void build_striplet_(Node *root);
        std::unordered_set<index_t> wg_edges(Node *root, Taxa &subset, weight_t ***graph);
        std::unordered_set<index_t> wb_edges(Node *root, Taxa &subset, weight_t ***graph);
        template <typename function1, typename function2>
        weight_t squartet(function1 f1, function2 f2, index_t size, index_t x, index_t y);
        void test_ssinglet(Node *root, Taxa &subset);
        void test_ssinglet_(Node *root, Taxa &subset);
        void test_sdoublet(Node *root, Taxa &subset);
        void test_sdoublet_(Node *root, Taxa &subset);
        void test_striplet(Node *root, Taxa &subset);
        void test_striplet_(Node *root, Taxa &subset);
        void test_pxlet(Node *root, std::unordered_set<index_t> &subtree, Taxa &subset);
        void test_pxlet_(Node *root, std::unordered_set<index_t> &subtree, Taxa &subset);
        void test_graph(Node *root, Taxa &subset, weight_t ***graph);
        void build_wstates(Node *root);
        void build_ssinglet(Node *root, std::unordered_map<index_t, index_t> quad);
        weight_t get_freq(std::unordered_map<index_t, index_t> quad);
        weight_t freq_(Node *root);
        void clear_wstates_(Node *root);
};

class SpeciesTree : public Tree {
    public:
        SpeciesTree(std::vector<Tree *> &input, Dict *dict, std::string mode, unsigned long int iter_limit);
        ~SpeciesTree();
        std::string annotate(std::vector<Tree *> input);
    private:
        index_t artifinyms;
        std::string mode;
        unsigned long int iter_limit;
        index_t artifinym();
        Node *construct_stree(std::vector<Tree *> &input, Taxa &subset, index_t parent_pid, index_t depth);
        Node *construct_stree(std::unordered_map<quartet_t, weight_t> &input, Taxa &subset, index_t parent_pid, index_t depth);
        Node *reroot(Node *root, std::unordered_set<index_t> &visited);
        Node *reroot_stree(Node *root, index_t artificial);
        Node *artificial2node(Node *root, index_t artificial);
        void get_freq(Node *root, std::vector<Tree *> input);
        std::string display_tree_annotated(Node *root);
};

extern std::ofstream subproblem_csv, quartets_txt, good_edges_txt, bad_edges_txt;
extern std::string verbose;
extern unsigned long long count[8];

#endif
