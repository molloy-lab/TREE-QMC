#include "tree.hpp"
#include "graph.hpp"

SpeciesTree::SpeciesTree(std::vector<Tree *> &input, Dict *dict, std::string mode, unsigned long int iter_limit, SpeciesTree *rst) {
    this->dict = dict;
    this->artifinyms = dict->max_size();
    this->mode = mode;
    this->iter_limit = iter_limit;
    this->root = NULL;
    Taxa subset(dict, mode);
    std::vector<Taxa> subsets;
    std::vector<Node *> internal = rst->decompose(subset, subsets);
    set_pseudonym(rst->pseudonyms);
    std::vector<Node *> subtrees;
    std::unordered_map<index_t, Node *> index2tree;
    for (index_t i = 0; i < subsets.size(); i ++) {
        Node *subtree = construct_stree(input, subsets[i], -1, 0);
        subtrees.push_back(subtree);
        index2tree[internal[i]->index] = subtree;
    }
    root = subtrees[0];
    std::unordered_set<Node *> visited;
    visited.insert(root);
    compose(root, subtrees, index2tree, visited, internal[0]->index);
    avail_pseudonyms.push_back(internal[0]->index);
    // std::cout << pseudonyms << " " << avail_pseudonyms.size() << std::endl;
}

void SpeciesTree::compose(Node *root, std::vector<Node *> &subtrees, std::unordered_map<index_t, Node *> &index2tree, std::unordered_set<Node *> &visited, index_t parent_index) {
    std::vector<Node *> leaves;
    get_leaves(root, &leaves);
    for (Node *leaf : leaves) {
        auto iter = index2tree.find(leaf->index);
        if (iter != index2tree.end()) {
            Node *subtree = iter->second;
            if (visited.find(subtree) == visited.end()) {
                visited.insert(subtree);
                compose(subtree, subtrees, index2tree, visited, leaf->index);
                Node *new_leaf = reroot_stree(subtree, parent_index);
                for (index_t i = 0; i < leaf->parent->children.size(); i ++) {
                    if (leaf->parent->children[i] == leaf) 
                        leaf->parent->children[i] = new_leaf;
                }
                new_leaf->parent = leaf->parent;
                avail_pseudonyms.push_back(leaf->index);
                delete leaf;
            }
        }
    }
}

std::vector<Node *> SpeciesTree::decompose(Taxa &subset, std::vector<Taxa> &subsets) {
    std::unordered_set<Node *> visited;
    traverse(root, visited, NULL, true);
    std::vector<Node *> internal;
    for (Node *node : visited) {
        if (node->children.size() != 0)
            internal.push_back(node);
    }
    for (Node *node : internal) {
        visited.clear();
        visited.insert(node);
        Taxa new_subset(subset);
        traverse(node, visited, &new_subset, true);
        subsets.push_back(new_subset);
    }
    return internal;
}

void SpeciesTree::traverse(Node *root, std::unordered_set<Node *> &visited, Taxa *subset, bool is_root) {
    std::vector<index_t> true_children;
    for (Node *child : root->children) {
        if (visited.find(child) == visited.end()) {
            visited.insert(child);
            traverse(child, visited, subset, false);
            true_children.push_back(child->index);
        }
    }
    if (root->parent != NULL && visited.find(root->parent) == visited.end()) {
        visited.insert(root->parent);
        traverse(root->parent, visited, subset, false);
        true_children.push_back(root->parent->index);
    }
    if (true_children.size() != 0) {
        if (! is_root && subset != NULL) 
            subset->struct_update(true_children, root->index);
    }
}
