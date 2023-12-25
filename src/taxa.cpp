#include "taxa.hpp"
#include "dict.hpp"

Taxa::Node::Node(index_t index) {
    parent = NULL;
    this->index = index;
    r_index = -1;
}

bool Taxa::Node::is_singleton() {
    return singleton;
}

Taxa::Taxa() {

}

Taxa::Taxa(Dict *dict, std::string mode) {
    this->updated = false;
    this->dict = dict;
    this->mode = mode;
    this->normal = mode[0];
    this->shared = mode[2];
    index2node = new Node*[dict->max_size()];
    std::memset(index2node, 0, sizeof(Node *) * dict->max_size());
    for (index_t i = 0; i < dict->size(); i ++) {
        Node *node = new Node(i);
        index2node[i] = node;
        node->singleton = true;
        roots.push_back(node);
        leaves.push_back(node);
    }
    this->singletons = roots.size();
}

Taxa::Taxa(const Taxa &taxa) {
    this->updated = false;
    singletons = taxa.singletons;
    mode = taxa.mode;
    normal = taxa.normal;
    shared = taxa.shared;
    dict = taxa.dict;
    index2node = new Node*[dict->max_size()];
    std::memset(index2node, 0, sizeof(Node *) * dict->max_size());
    for (index_t i = 0; i < dict->max_size(); i ++) {
        if (taxa.index2node[i] == NULL) continue;
        Node *new_node = new Node(i);
        index2node[i] = new_node;
        new_node->r_index = taxa.index2node[i]->r_index;
        new_node->singleton = taxa.index2node[i]->singleton;
    }
    for (index_t i = 0; i < dict->max_size(); i ++) {
        if (index2node[i] == NULL) continue;
        Node *new_node = index2node[i];
        if (taxa.index2node[i]->parent == NULL) 
            new_node->parent = NULL;
        else 
            new_node->parent = index2node[taxa.index2node[i]->parent->index];
    }
    for (Node *root : taxa.roots) 
        roots.push_back(index2node[root->index]);
    for (Node *leaf : taxa.leaves) 
        leaves.push_back(index2node[leaf->index]);
}

Taxa::~Taxa() {
    for (index_t i = 0; i < dict->max_size(); i ++) {
        if (index2node[i] != NULL) delete index2node[i];
    }
    delete [] index2node;
}

void Taxa::struct_update(std::vector<index_t> &subset, index_t artificial) {
    Node *new_root = new Node(artificial);
    index2node[artificial] = new_root;
    roots.push_back(new_root);
    new_root->singleton = false;
    for (index_t index : subset) {
        Node *node = index2node[index];
        node->parent = new_root;
        node->r_index = -1;
        node->singleton = false;
    }
    sort_taxa();
}

void Taxa::weight_update(std::unordered_map<index_t, index_t> &subset) {
    if (shared == '1') {
        if (! updated) {
            updated = true;
            // std::cout << size() << std::endl;
            for (Node *node : leaves) {
                node->singleton = node->parent == NULL;
            }
            if (normal == '1' || normal == '2') {
                std::queue<Node *> queue;
                std::unordered_set<Node *> visited;
                for (Node *node : leaves) {
                    queue.push(node);
                    visited.insert(node);
                    node->degree = 1;
                }
                while (! queue.empty()) {
                    Node *head = queue.front();
                    queue.pop();
                    if (head->parent != NULL) {
                        if (visited.find(head->parent) == visited.end()) {
                            queue.push(head->parent);
                            visited.insert(head->parent);
                            head->parent->degree = 0;
                            head->parent->size = 0;
                        }
                        head->parent->degree += 1;
                    }
                }
                if (normal == '1') {
                    for (Node *node : leaves) {
                        queue.push(node);
                        node->size = 1.0;
                    }
                    while (! queue.empty()) {
                        Node *head = queue.front();
                        queue.pop();
                        if (head->parent != NULL) {
                            head->parent->degree -= 1;
                            head->parent->size += head->size;
                            if (head->parent->degree == 0) {
                                queue.push(head->parent);
                            }
                        }
                    }
                }
            }
            sort_taxa();
        }
    }
    else {
        for (Node *node : leaves) {
            if (subset.find(node->index) != subset.end()) {
                node->singleton = subset[node->index] == 1 && node->parent == NULL;
            }
        }
        if (normal == '1' || normal == '2') {
            std::queue<Node *> queue;
            std::unordered_set<Node *> visited;
            for (Node *node : leaves) {
                if (subset.find(node->index) != subset.end()) {
                    queue.push(node);
                    visited.insert(node);
                    node->degree = subset[node->index];
                }
            }
            while (! queue.empty()) {
                Node *head = queue.front();
                queue.pop();
                if (head->parent != NULL) {
                    if (visited.find(head->parent) == visited.end()) {
                        queue.push(head->parent);
                        visited.insert(head->parent);
                        head->parent->degree = 0;
                        head->parent->size = 0;
                    }
                    head->parent->degree += 1;
                }
            }
            if (normal == '1') {
                for (Node *node : leaves) {
                    if (subset.find(node->index) != subset.end()) {
                        queue.push(node);
                        node->size = subset[node->index];
                    }
                }
                while (! queue.empty()) {
                    Node *head = queue.front();
                    queue.pop();
                    if (head->parent != NULL) {
                        head->parent->degree -= 1;
                        head->parent->size += head->size;
                        if (head->parent->degree == 0) {
                            queue.push(head->parent);
                        }
                    }
                }
            }
        }
        sort_taxa();
    }
}

std::string Taxa::to_list() {
    std::string s = "";
    std::vector<std::string> root_indices;
    for (Node *root : roots) {
        root_indices.push_back(dict->index2label(root->index));
    }
    std::sort(root_indices.begin(), root_indices.end());
    for (std::string index : root_indices) {
        s += index + " ";
    }
    return s;
}

std::string Taxa::to_string() {
    std::string s = "";
    for (Node *node : leaves) {
        s += "(" + std::to_string((double) root_weight(node->index)) + ")";
        while (node != NULL) {
            s += std::to_string(node->index) + ":" + std::to_string(node->degree) + "," + std::to_string(node->size) + "," + std::to_string(node->singleton) + " ";
            node = node->parent;
        }
        s += "\n";
    }
    for (Node *root : roots) {
        s += dict->index2label(root->index) + "/" + std::to_string(root->index) + ":" + std::to_string(root->degree) + "," + std::to_string(root->size) + "," + std::to_string(root->singleton) + " ";
    }
    return s + "\n";
}

index_t Taxa::size() {
    return roots.size();
}

char Taxa::normalization() {
    return normal;
}

char Taxa::get_shared() {
    return shared;
}

weight_t Taxa::get_sum() {
    weight_t count = roots.size() - 2;
    return (count - 1) * count;
}

bool Taxa::is_singleton(index_t index) {
    return index2node[index]->is_singleton();
}

index_t Taxa::singleton_taxa() {
    return singletons;
}

index_t Taxa::artificial_taxa() {
    return roots.size() - singletons;
}

index_t Taxa::leaf_at(index_t i) {
    return leaves[i]->index;
}

index_t Taxa::root_at(index_t i) {
    return roots[i]->index;
}

index_t Taxa::artificial_at(index_t i) {
    return roots[i - 1 + singletons]->index;
}

index_t Taxa::get_index(index_t index) {
    return get_root(index)->index;
}

index_t Taxa::root_index(index_t index) {
    return get_root(index)->r_index;
}

index_t Taxa::root_key(index_t index) {
    Node *root = get_root(index);
    if (root->is_singleton()) return 0;
    return root->r_index - singletons + 1;
}

Taxa::Node *Taxa::get_root(index_t index) {
    Node *root = index2node[index];
    return get_root(root);
}

Taxa::Node *Taxa::get_root(Taxa::Node *root) {
    while (root->parent != NULL) 
        root = root->parent;
    return root;
}

weight_t Taxa::root_weight(index_t index) {
    if (normal == '0') {
        return 1.0;
    }
    else if (normal == '1') {
        return 1.0 / get_root(index)->size;
    }
    else {
        return get_weight(index2node[index]);
    }
}

weight_t Taxa::get_weight(Node *root) {
    weight_t w = 1.0 / root->degree;
    while (root->parent != NULL) {
        w /= (weight_t) root->parent->degree;
        root = root->parent;
    }
    return w;
}

void Taxa::sort_taxa() {
    std::vector<Node *> new_roots = std::vector<Node *>();
    for (Node *root : roots) {
        if (root->is_singleton()) {
            new_roots.push_back(root);
        }
    }
    singletons = new_roots.size();
    for (Node *root : roots) {
        if (! root->is_singleton() && root->parent == NULL) {
            new_roots.push_back(root);
        }
    }
    roots.clear();
    index_t index = 0;
    for (Node *root : new_roots) {
        roots.push_back(root);
        root->r_index = index ++;
    }
}
