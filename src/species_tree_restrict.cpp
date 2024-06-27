#include "tree.hpp"
#include "graph.hpp"

SpeciesTree::SpeciesTree(std::vector<std::vector<index_t>> &clades, std::vector<std::string> &names, Dict *dict) {
    this->dict = dict;
    this->root = NULL;
    index_t row = dict->size(), col = clades.size();
    index_t **matrix = new index_t*[row];
    for (index_t i = 0; i < row; i ++) {
        matrix[i] = new index_t[col];
        for (index_t j = 0; j < col; j ++) 
            matrix[i][j] = 0;
    }
    for (index_t j = 0; j < col; j ++) {
        for (index_t i : clades[j]) 
            matrix[i][j] = 1;
    }
    for (index_t i = 0; i < col; i ++) {
        for (index_t j = i + 1; j < col; j ++) {
            index_t k = 0;
            for (k = 0; matrix[k][i] == matrix[k][j]; k ++) ;
            if (matrix[k][i] < matrix[k][j]) {
                for (k = 0; k < row; k ++) {
                    index_t temp = matrix[k][i];
                    matrix[k][i] = matrix[k][j];
                    matrix[k][j] = temp;
                }
                std::string temps = names[i];
                names[i] = names[j];
                names[j] = temps;
                std::vector<index_t> tempc = clades[i];
                clades[i] = clades[j];
                clades[j] = tempc;
            }
        }
    }
    index_t **lower = new index_t*[row];
    for (index_t i = 0; i < row; i ++) {
        lower[i] = new index_t[col];
        for (index_t j = 0; j < col; j ++) 
            lower[i][j] = -2;
    }
    for (index_t i = 0; i < row; i ++) {
        index_t prev = -1;
        for (index_t j = 0; j < col; j ++) {
            if (matrix[i][j] == 1) {
                lower[i][j] = prev;
                prev = j;
            }
        }
    }
    bool compatible = true;
    index_t *max = new index_t[col];
    for (index_t j = 0; j < col; j ++) {
        max[j] = -2;
        for (index_t i = 0; i < row; i ++) {
            if (max[j] < lower[i][j]) 
                max[j] = lower[i][j];
        }
        //std::cout << max[j] << " ";
        for (index_t i = 0; i < row; i ++) {
            if (matrix[i][j] == 1 && lower[i][j] != max[j]) 
                compatible = false;
        }
    }
    if (! compatible) {
        std::cout << "\nWARNING: Incompatible clades" << std::endl;
        exit(1);
    }
    else {
        std::cout << "Input clades are compatible" << std::endl;
        std::cout << "Building perfect phylogeny" << std::endl;
    }
    //std::unordered_map<index_t, std::string> index2clade;
    root = new Node(pseudonym());
    std::vector<Node *> internal;
    for (index_t j = 0; j < col; j ++) {
        Node *node = new Node(pseudonym());
        internal.push_back(node);
        //index2clade[node->index] = names[j];
        //std::cout << names[j] << " " << node->index << std::endl;
    }
    for (index_t j = 0; j < col; j ++) {
        Node *parent = max[j] == -1 ? root : internal[max[j]];
        internal[j]->parent = parent;
        parent->children.push_back(internal[j]);
    }
    for (index_t i = 0; i < row; i ++) {
        index_t j;
        for (j = col - 1; matrix[i][j] == 0; j --) ;
        Node *node = internal[j];
        Node *new_node = new Node(i);
        node->children.push_back(new_node);
        new_node->parent = node;
    }
    root = simplify(root);
    for (index_t i = 0; i < row; i ++) 
        delete [] matrix[i];
    delete [] matrix;
    for (index_t i = 0; i < row; i ++) 
        delete [] lower[i];
    delete [] lower;
}

Node *SpeciesTree::simplify(Node *root) {
    if (root->children.size() == 1) {
        Node *child = root->children[0];
        root->children.clear();
        avail_pseudonyms.push_back(root->index);
        delete root;
        child->parent = NULL;
        return simplify(child);
    }
    else {
        for (index_t i = 0; i < root->children.size(); i ++) {
            root->children[i] = simplify(root->children[i]);
            root->children[i]->parent = root;
        }
        return root;
    }
}

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
    std::vector<Node *> internal;
    get_internal_nodes(root, &internal);
    for (Node *node : internal) {
        std::unordered_set<Node *> visited;
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
