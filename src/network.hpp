#if ENABLE_TOB
#ifndef __NETWORK__
#define __NETWORK__

#include "tree.hpp"
#include "rlib_dirs.hpp"





class NetworkNode {
    friend class Network;
    public: 
        NetworkNode(std::string label);
        ~NetworkNode();
        static std::size_t get_count();
        double get_lambda() {
            return lambda;
        }
        std::size_t get_degree() {
            return children.size() + 1;
        }
    private: 
        NetworkNode *parent, *hybrid;
        std::vector<NetworkNode*> children;
        std::string label;
        double length, support, lambda;
        bool bridge;
        static std::size_t count;
};

class Network {
    public: 
        Network(std::string &newick);
        Network(Network &network);
        Network(Tree * hybrid_tree, Dict *dict);
        // Network(Node *root, Dict *dict);
        ~Network();
        std::string to_string();
        std::string to_string_basic();
        void compress(double coeffi);
        std::vector<NetworkNode *> get_hybridization();
        std::vector<NetworkNode *> get_blobs();
        std::unordered_map<std::string, NetworkNode*> get_label2hybrid();
        Node* reroot_hybrid_tree_on_edge_above_node(Node* node, Tree* tre);
        std::size_t get_nodes() {
            std::vector<NetworkNode *> nodes = get_nodes(root);
            std::size_t count = 0;
            for (auto node : nodes) {
                if (node->children.size() == 0) continue;
                if (node->parent == NULL) continue;
                if (node->parent->parent == NULL && node->parent->children.size() == 2 && node != node->parent->children[1]) continue;
                count ++;
            }
            return count;
        }
        
    private: 
        size_t pseudonyms;
        static size_t hybrid_pseudonyms;
        NetworkNode *root;
        Tree *tob_annotation_tree;
        std::unordered_map<std::string, NetworkNode*> label2node, label2hybrid;
        static std::string gen_hybrid_id();
        NetworkNode *build_network(Node * root, Dict *dict, std::unordered_map<std::string, NetworkNode*> &label2hybrid);
        NetworkNode *build_network(std::string newick, std::unordered_map<std::string, NetworkNode*> &label2hybrid);
        bool check_root(Node *root, std::unordered_set<index_t> &candidates);
        std::string display_network(NetworkNode *root);
        std::string display_network_basic(NetworkNode *root);
        std::string display_network_bridge(NetworkNode *root);
        std::string pseudonym();
        void depth_first_search(NetworkNode *root, NetworkNode *parent, std::size_t depth, std::unordered_map<NetworkNode *, std::pair<std::size_t, std::size_t>> &visited);
        NetworkNode *build_tree_of_blob(NetworkNode *root);
        NetworkNode *simplify_tree(NetworkNode *root);
        void compress_length(NetworkNode *root, double coeffi);
        std::vector<NetworkNode *> get_hybridization(NetworkNode *root);
        std::vector<NetworkNode *> get_blobs(NetworkNode *root);
        size_t find_child_index_in_children(Node* node, Node * child2find);
        void update_hybird_info(Node* node, index_t old_child_index);
        Node* rerooting_hybrid_tree_above_node(Node * node, Node * child);
        std::pair<size_t, std::vector<size_t>> find_above_index_in_cycle(Node *root);
        std::vector<NetworkNode *> get_nodes(NetworkNode *root) {
            std::vector<NetworkNode *> nodes;
            nodes.push_back(root);
            for (auto child : root->children) {
                std::vector<NetworkNode *> sub_nodes = get_nodes(child);
                for (auto node : sub_nodes) {
                    nodes.push_back(node);
                }
            }
            return nodes;
        }
        

};

#endif // __NETWORK__  // ENABLE_TOB
#endif