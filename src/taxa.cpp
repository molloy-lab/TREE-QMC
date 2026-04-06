#include "taxa.hpp"
#include "dict.hpp"

Taxa::Node::Node(index_t index) {
    parent = NULL;
    this->index = index;
    r_index = -1;
    visit_token = 0;
}

bool Taxa::Node::is_singleton() {
    return singleton;
}

Taxa::Taxa() {
    updated = false;
    cache_valid = false;
    visit_epoch = 0;
}

Taxa::Taxa(Dict *dict, std::string mode) {
    this->updated = false;
    this->cache_valid = false;
    this->visit_epoch = 0;
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
    cached_root_index.assign(dict->max_size(), -1);
    cached_root_rindex.assign(dict->max_size(), -1);
    cached_root_key.assign(dict->max_size(), -1);
    cached_root_weight.assign(dict->max_size(), 0.0);
    traversal_queue.reserve(dict->max_size());
    this->singletons = roots.size();
}

Taxa::Taxa(const Taxa &taxa) {
    this->updated = false;
    this->cache_valid = false;
    this->visit_epoch = 0;
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
    cached_root_index.assign(dict->max_size(), -1);
    cached_root_rindex.assign(dict->max_size(), -1);
    cached_root_key.assign(dict->max_size(), -1);
    cached_root_weight.assign(dict->max_size(), 0.0);
    traversal_queue.reserve(dict->max_size());
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
    cache_valid = false;
    sort_taxa();
}

void Taxa::rebuild_lookup_cache() {
    const index_t max_nodes = dict->max_size();
    std::fill(cached_root_index.begin(), cached_root_index.end(), static_cast<index_t>(-1));
    std::fill(cached_root_rindex.begin(), cached_root_rindex.end(), static_cast<index_t>(-1));
    std::fill(cached_root_key.begin(), cached_root_key.end(), static_cast<index_t>(-1));
    std::fill(cached_root_weight.begin(), cached_root_weight.end(), 0.0);
    for (index_t i = 0; i < max_nodes; i ++) {
        Node *node = index2node[i];
        if (node == NULL) continue;
        Node *root = get_root(node);
        cached_root_index[i] = root->index;
        cached_root_rindex[i] = root->r_index;
        cached_root_key[i] = root->is_singleton() ? 0 : root->r_index - singletons + 1;
        if (normal == '0') {
            cached_root_weight[i] = 1.0;
        }
        else if (normal == '1') {
            cached_root_weight[i] = 1.0 / root->size;
        }
        else {
            cached_root_weight[i] = get_weight(node);
        }
    }
    cache_valid = true;
}

void Taxa::weight_update(std::unordered_map<index_t, index_t> &subset) {
    if (shared == '1') {
        if (! updated) {
            updated = true;
            for (Node *node : leaves) {
                node->singleton = node->parent == NULL;
            }
            if (normal == '1' || normal == '2') {
                traversal_queue.clear();
                if (visit_epoch == std::numeric_limits<unsigned long long>::max()) {
                    visit_epoch = 0;
                    for (index_t i = 0; i < dict->max_size(); i ++) {
                        if (index2node[i] != NULL) index2node[i]->visit_token = 0;
                    }
                }
                const unsigned long long token = ++visit_epoch;
                for (Node *node : leaves) {
                    traversal_queue.push_back(node);
                    node->visit_token = token;
                    node->degree = 1;
                }
                for (std::size_t qhead = 0; qhead < traversal_queue.size(); ++qhead) {
                    Node *head = traversal_queue[qhead];
                    if (head->parent != NULL) {
                        if (head->parent->visit_token != token) {
                            traversal_queue.push_back(head->parent);
                            head->parent->visit_token = token;
                            head->parent->degree = 0;
                            head->parent->size = 0;
                        }
                        head->parent->degree += 1;
                    }
                }
                if (normal == '1') {
                    traversal_queue.clear();
                    for (Node *node : leaves) {
                        traversal_queue.push_back(node);
                        node->size = 1.0;
                    }
                    for (std::size_t qhead = 0; qhead < traversal_queue.size(); ++qhead) {
                        Node *head = traversal_queue[qhead];
                        if (head->parent != NULL) {
                            head->parent->degree -= 1;
                            head->parent->size += head->size;
                            if (head->parent->degree == 0) {
                                traversal_queue.push_back(head->parent);
                            }
                        }
                    }
                }
            }
            sort_taxa();
            rebuild_lookup_cache();
        }
    }
    else {
        for (Node *node : leaves) {
            auto it = subset.find(node->index);
            if (it != subset.end()) {
                node->singleton = it->second == 1 && node->parent == NULL;
            }
        }
        if (normal == '1' || normal == '2') {
            traversal_queue.clear();
            if (visit_epoch == std::numeric_limits<unsigned long long>::max()) {
                visit_epoch = 0;
                for (index_t i = 0; i < dict->max_size(); i ++) {
                    if (index2node[i] != NULL) index2node[i]->visit_token = 0;
                }
            }
            const unsigned long long token = ++visit_epoch;
            for (Node *node : leaves) {
                auto it = subset.find(node->index);
                if (it != subset.end()) {
                    traversal_queue.push_back(node);
                    node->visit_token = token;
                    node->degree = it->second;
                }
            }
            for (std::size_t qhead = 0; qhead < traversal_queue.size(); ++qhead) {
                Node *head = traversal_queue[qhead];
                if (head->parent != NULL) {
                    if (head->parent->visit_token != token) {
                        traversal_queue.push_back(head->parent);
                        head->parent->visit_token = token;
                        head->parent->degree = 0;
                        head->parent->size = 0;
                    }
                    head->parent->degree += 1;
                }
            }
            if (normal == '1') {
                traversal_queue.clear();
                for (Node *node : leaves) {
                    auto it = subset.find(node->index);
                    if (it != subset.end()) {
                        traversal_queue.push_back(node);
                        node->size = it->second;
                    }
                }
                for (std::size_t qhead = 0; qhead < traversal_queue.size(); ++qhead) {
                    Node *head = traversal_queue[qhead];
                    if (head->parent != NULL) {
                        head->parent->degree -= 1;
                        head->parent->size += head->size;
                        if (head->parent->degree == 0) {
                            traversal_queue.push_back(head->parent);
                        }
                    }
                }
            }
        }
        sort_taxa();
        rebuild_lookup_cache();
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
    if (cache_valid) return cached_root_index[index];
    return get_root(index)->index;
}

index_t Taxa::root_index(index_t index) {
    if (cache_valid) return cached_root_rindex[index];
    return get_root(index)->r_index;
}

index_t Taxa::root_key(index_t index) {
    if (cache_valid) return cached_root_key[index];
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
    if (cache_valid) return cached_root_weight[index];
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
