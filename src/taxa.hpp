#ifndef TAXA_HPP
#define TAXA_HPP

#include "utility.hpp"
#include "dict.hpp"

class Taxa {
    public:
        class Node {
            friend class Taxa;
            public:
                Node(index_t index);
                bool is_singleton();
            private:
                Node *parent;
                index_t index, r_index, size, degree;
                bool singleton;
        };
        Taxa();
        Taxa(Dict *dict, std::string mode);
        Taxa(const Taxa &taxa);
        ~Taxa();
        void struct_update(std::vector<index_t> &subset, index_t artificial);
        void weight_update(std::unordered_map<index_t, index_t> &subset);
        std::string to_list();
        std::string to_string();
        index_t size();
        char normalization();
        char get_shared();
        weight_t get_sum();
        index_t singleton_taxa();
        index_t artificial_taxa();
        bool is_singleton(index_t index);
        index_t leaf_at(index_t i);
        index_t root_at(index_t i);
        index_t artificial_at(index_t i);
        index_t get_index(index_t index);
        index_t root_index(index_t index);
        index_t root_key(index_t index);
        weight_t root_weight(index_t index);
        Node *get_root(index_t index);
        Node *get_root(Node *root);
        weight_t get_weight(Node *root);
    private:
        std::vector<Node *> leaves, roots;
        // Node **index2node;
        std::unordered_map<index_t, Node*> index2node;
        Dict *dict;
        index_t singletons;
        std::string mode;
        char normal, shared;
        bool updated;
        void sort_taxa();
};

#endif
