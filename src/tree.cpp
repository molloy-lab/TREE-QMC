#include "tree.hpp"

Tree::Tree() {
    pseudonyms = 0;
}

Tree::Tree(const std::string &newick,
           Dict *dict,
           const std::unordered_map<std::string, std::string> &indiv2taxon,
           weight_t support_default) {
    pseudonyms = 0;
    this->dict = dict;
    this->support_default = support_default;
    this->pcs_node = NULL;
    root = build_tree(newick, indiv2taxon);
}

Tree::~Tree() {
    // TODO: Do we need to check this?
    delete root;
}

Node* Tree::get_root() {
    return root;
}

std::string Tree::to_string() {
    return display_tree(root) + ";";
}

std::string Tree::to_string_basic() {
    return display_tree_basic(root) + ";";
}

size_t Tree::refine() {
    bool flag = root->children.size() == 3;
    size_t count = refine_tree(root);
    if (flag) count -= 1;
    return count;
}

void Tree::prepare(std::string weight_mode, weight_t low, weight_t high, weight_t threshold) {
    //std::cout << display_tree(root) << std::endl;
    prepare_tree(root, weight_mode, low, high, threshold);
}

index_t Tree::size() {
    return index2node.size();
}

std::unordered_map<index_t, index_t> &Tree::get_indices() {
    return indices;
}

weight_t ***Tree::build_graph(Taxa &subset) {
    index_t s = subset.singleton_taxa(), m = subset.artificial_taxa();
    build_states(root, subset);
    for (index_t i = 0; i <= m; i ++) 
        sa_doublet(root, 0, i);
    for (index_t i = 1; i <= m; i ++) {
        for (index_t j = 0; j <= m; j ++) {
            if (i == j) continue;
            aa_doublet(root, i, j);
        }
    }
    weight_t ***graph = new weight_t**[2];
    graph[0] = Matrix::new_mat(subset.size());
    graph[1] = Matrix::new_mat(subset.size());
    bad_edges(root, subset, graph);
    good_edges(root, subset, graph);
    clear_states(root);
    return graph;
}

index_t Tree::pseudonym() {
    return - (++ pseudonyms);
}

Node *Tree::build_tree(const std::string &newick,
                       const std::unordered_map<std::string, std::string> &indiv2taxon) {
    std::string label, leaf_label, support, length;
    std::size_t sep;

    if (newick.length() == 0 || newick.at(0) != '(') {
        sep = newick.find(":");
        if (sep != std::string::npos) {
            leaf_label = newick.substr(0, sep);
            length = newick.substr(sep + 1, std::string::npos);
            //std::cout << "Found leaf node " << leaf_label << ' ' << length << std::endl;
        } else {
            leaf_label = newick.substr(0, std::string::npos);
            //std::cout << "Found leaf node " << leaf_label << std::endl;
        }

        // Map leaf label if possible
        auto itr = indiv2taxon.find(leaf_label);
        if (itr == indiv2taxon.end())
            label = leaf_label;
        else
            label = itr->second;

        // Create leaf node
        Node *root = new Node(dict->label2index(label));
        if (length.size() > 0) root->length = std::stod(length);

        // Update related data structures
        index2node[root->index] = root;
        if (indices.find(root->index) == indices.end())
            indices[root->index] = 0;
        indices[root->index] ++;

        return root;
    }
    else {
        // Create internal node
        Node *root = new Node(pseudonym());

        // Recurse on all but last child
        int k = 1;
        for (int i = 0, j = 0; i < newick.length(); i ++) {
            if (newick.at(i) == '(') j ++;
            if (newick.at(i) == ')') j --;
            if (newick.at(i) == ',' && j == 1) {
                root->children.push_back(build_tree(newick.substr(k, i - k), indiv2taxon));
                k = i + 1;
            }
        }

        int i = newick.length() - 1;
        while (newick.at(i) != ')') i --;
        std::string branch = newick.substr(i + 1, std::string::npos);
        if (branch.find(";") == std::string::npos) {
            // Get internal branch and support
            
            sep = branch.find(":");
            if (sep != std::string::npos) {
                support = branch.substr(0, sep);
                length = branch.substr(sep + 1, std::string::npos);
                //std::cout << "Found internal node " << branch << " with support " << support << " and length " << length << std::endl;
                if (length.size() > 0) root->length = std::stod(length);
            } else {
                support = branch.substr(0, std::string::npos);
                //std::cout << "Found internal node " << branch << ' ' << support << std::endl;
            }
            if (support.size() > 0) {
                if (support == "PCS")
                    pcs_node = root;
                else
                    root->support = std::stod(support);
            }
            else root->support = support_default;  // allows user to change default support
        }
        //else {
        //    std::cout << "Found root" << std::endl;
        //}

        // Recurse on last child
        root->children.push_back(build_tree(newick.substr(k, i - k), indiv2taxon));

        // Set root as children's parent
        for (Node *child : root->children)
            child->parent = root;

        return root;
    }
}


std::string Tree::display_tree(Node *root) {
    if (root->children.size() == 0) 
        return dict->index2label(root->index) + ":" + std::to_string((double) root->length);
    std::string s = "(";
    for (Node * node : root->children) 
        s += display_tree(node) + ",";
    s[s.size() - 1] = ')';
    return s + std::to_string((double) root->support) + ":" + std::to_string((double) root->length);
}

std::string Tree::display_tree_basic(Node *root) {
    if (root->children.size() == 0) 
        return dict->index2label(root->index);
    std::string s = "(";
    for (Node * node : root->children) 
        s += display_tree_basic(node) + ",";
    s[s.size() - 1] = ')';
    return s;
}

std::string Tree::display_tree_index(Node *root) {
    if (root->children.size() == 0) 
        return std::to_string(root->index);
    std::string s = "(";
    for (Node * node : root->children) 
        s += display_tree_index(node) + ",";
    s[s.size() - 1] = ')';
    return s;
}

size_t Tree::refine_tree(Node *root) {
    size_t total = 0;

    for (index_t i = 0; i < root->children.size(); i ++) 
        total += refine_tree(root->children[i]);

    while (root->children.size() > 2) {
        index_t i = rand() % root->children.size(), j = i;
        while (j == i) j = rand() % root->children.size();
        Node *new_root = new Node(pseudonym());
        new_root->support = 0.0;  //Node *new_root = new Node(pseudonym(), true);
        total ++;
        new_root->children.push_back(root->children[i]);
        new_root->children.push_back(root->children[j]);
        root->children[i]->parent = root->children[j]->parent = new_root;
        root->children.erase(root->children.begin() + i);
        root->children.erase(root->children.begin() + (j > i ? j - 1 : j));
        root->children.push_back(new_root);
        new_root->parent = root;
    }

    return total;
}

void Tree::prepare_tree(Node *root, std::string weight_mode, weight_t low, weight_t high, weight_t threshold) {
    //assert(root->children.size() == 0 || root->children.size() == 2);

    weight_t s = root->support;

    // Map support to 0 and 1 interval
    if (s < low || s > high) s = low;
        s = (s - low) / (high - low);

    // Contract low support branches by setting support to 0
    if (s < threshold) s = 0;

    // If not weighting by branch support, set support values to 1
    // do for fast mode as well because will be used for support estimation
    // this is already done in instance.cpp
    if (weight_mode == "n" || weight_mode == "f" || weight_mode == "l")
        if (s > 0) s = 1;

    // Save these values
    root->support = s;
    root->support_[0] = 1 - s;
    root->support_[1] = 1;

    // Handle length values
    if (weight_mode == "h" || weight_mode == "l")
        root->length_ = exp(- root->length);
    else
        root->length_ = 1;

    // Continue
    for (Node *child : root->children)
        prepare_tree(child, weight_mode, low, high, threshold);
}

void Tree::clear_states(Node *root) {
    root->delete_states();
    for (Node *child : root->children) 
        clear_states(child);
}

void Tree::build_states(Node *root, Taxa &subset) { 
    root->new_states(subset.artificial_taxa());
    if (root->children.size() == 0) {
        index_t index = subset.root_key(root->index);
        root->singlet[index] = subset.root_weight(root->index);
    }
    else {
        for (Node *child : root->children) {
            build_states(child, subset);
            for (index_t i = 0; i < root->size + 1; i ++) 
                root->singlet[i] += child->singlet[i];
        }
    }
    root->s1 = root->s2 = 0;
    for (int i = 1; i < root->size + 1; i ++) {
        root->s1 += root->singlet[i];
        root->s2 += root->singlet[i] * root->singlet[i];
    }
}

weight_t Tree::get_doublet(Node *subtree, index_t x, index_t y, bool complement) {
    if (complement) {
        weight_t *c = new weight_t[subtree->size + 1], s1 = 0, s2 = 0;
        for (int i = 0; i <= subtree->size; i ++) {
            c[i] = root->singlet[i] - subtree->singlet[i];
            if (i != 0) {
                s1 += c[i];
                s2 += c[i] * c[i];
            }
        }
        weight_t ret = Node::get_doublet(c, s1, s2, x, y);
        delete [] c;
        return ret;
    }
    return Node::get_doublet(subtree->singlet, subtree->s1, subtree->s2, x, y);
}

void Tree::sa_doublet(Node *root, weight_t sum, index_t x) {
    root->add_doublet(0, x, sum);
    // root->doublet[0][x] = sum;
    if (root->children.size() == 0) return ;
    sa_doublet(root->children[0], sum + get_doublet(root->children[1], x, 0, false), x);
    sa_doublet(root->children[1], sum + get_doublet(root->children[0], x, 0, false), x);
}

/*
weight_t Tree::aa_doublet(Node *root, index_t x, index_t y) {
    if (root->children.size() == 0) return 0;
    weight_t t0 = get_doublet(root->children[0], x, y, false), t1 = get_doublet(root->children[1], x, y, false);
    weight_t s0 = t0 == 0 ? 0 : aa_doublet(root->children[0], x, y);
    weight_t s1 = t1 == 0 ? 0 : aa_doublet(root->children[1], x, y);
    root->doublet[x][y] = s0 + s1 + root->children[0]->singlet[x] * t1 + root->children[1]->singlet[x] * t0;
    if (root->doublet[x][y] == 0) total_count[3] ++;
    return root->doublet[x][y];
}
*/

weight_t Tree::aa_doublet(Node *root, index_t x, index_t y) {
    if (root->children.size() == 0) return 0;
    weight_t c0 = root->children[0]->singlet[x], c1 = root->children[1]->singlet[x], s0 = 0, s1 = 0;
    weight_t t0 = get_doublet(root->children[0], x, y, false), t1 = get_doublet(root->children[1], x, y, false);
    if (verbose > "1" || (c0 != 0 && t0 != 0)) s0 = aa_doublet(root->children[0], x, y);
    if (verbose > "1" || (c1 != 0 && t1 != 0)) s1 = aa_doublet(root->children[1], x, y);
    weight_t c = s0 + s1 + c0 * t1 + c1 * t0;
    if (verbose > "1") {
        count[1] ++;
        if (root->singlet[x] == 0 || get_doublet(root, x, y, false) == 0) count[2] ++;
        if (c < 1e-9) count[3] ++;
    }
    root->add_doublet(x, y, c);
    return c;
}

std::unordered_set<index_t> Tree::bad_edges(Node *root, Taxa &subset, weight_t ***graph) {
    if (root->children.size() == 0) {
        std::unordered_set<index_t> subtree;
        index_t index = subset.get_index(root->index);
        subtree.insert(index);
        return subtree;
    }
    else {
        std::unordered_set<index_t> l = bad_edges(root->children[0], subset, graph);
        std::unordered_set<index_t> r = bad_edges(root->children[1], subset, graph);
        std::unordered_set<index_t> subtree;
        for (auto i = l.begin(); i != l.end(); i ++) {
            if (subtree.find(*i) == subtree.end())
                subtree.insert(*i);
        }
        for (auto j = r.begin(); j != r.end(); j ++) {
            if (subtree.find(*j) == subtree.end()) 
                subtree.insert(*j);
        }
        for (auto i = l.begin(); i != l.end(); i ++) {
            for (auto j = r.begin(); j != r.end(); j ++) {
                if (*i == *j) continue;
                index_t x = subset.root_key(*i), y = subset.root_key(*j);
                index_t i_ = subset.root_index(*i), j_ = subset.root_index(*j);
                weight_t s = 0;
                if (x == 0) {
                    if (y == 0) {
                        s += index2node[*i]->get_doublet(0, 0) + index2node[*j]->get_doublet(0, 0);
                        s -= root->children[0]->get_doublet(0, 0) + root->children[1]->get_doublet(0, 0);
                        s += get_doublet(root, 0, 0, true);
                    }
                    else {
                        s += root->children[1]->get_doublet(y, 0);
                        weight_t t = root->children[1]->singlet[y];
                        s += t * (index2node[*i]->get_doublet(0, y) - root->children[0]->get_doublet(0, y));
                        s += t * get_doublet(root, y, 0, true);
                    }
                }
                else {
                    if (y == 0) {
                        s += root->children[0]->get_doublet(x, 0);
                        weight_t t = root->children[0]->singlet[x];
                        s += t * (index2node[*j]->get_doublet(0, x) - root->children[1]->get_doublet(0, x));
                        s += t * get_doublet(root, x, 0, true);
                    }
                    else {
                        weight_t t0 = root->children[0]->singlet[x], t1 = root->children[1]->singlet[y];
                        s += t1 * root->children[0]->get_doublet(x, y) + t0 * root->children[1]->get_doublet(y, x);
                        s += t1 * t0 * get_doublet(root, x, y, true);
                    }
                }
                graph[1][i_][j_] += s;
                graph[1][j_][i_] += s;
            }
        }
        return subtree;
    }
}

void Tree::good_edges(Node *root, Taxa &subset, weight_t ***graph) {
    index_t s = subset.singleton_taxa(), m = subset.artificial_taxa();
    std::unordered_map<index_t, index_t> tree_indices = get_indices();
    std::unordered_set<index_t> valid;
    for (auto elem : tree_indices) 
        valid.insert(subset.root_index(elem.first));
    if (valid.size() < 4) return ;
    if (subset.normalization() != '0') {
        weight_t *c = new weight_t[m + 1];
        for (index_t i = 0; i <= m; i ++) 
            c[i] = root->singlet[i];
        c[0] -= 2;
        weight_t sum = Node::get_doublet(c, root->s1, root->s2, 0, 0);
        delete [] c;
        for (index_t i = 0; i < subset.size(); i ++) {
            if (valid.find(i) == valid.end()) continue;
            for (index_t j = 0; j < subset.size(); j ++) {
                if (i == j || valid.find(j) == valid.end()) continue;
                graph[0][i][j] = sum - graph[1][i][j];
            }
        }
    }
    else {
        weight_t *c = new weight_t[m + 1];
        for (index_t i = 0; i <= m; i ++) 
            c[i] = root->singlet[i];
        c[0] -= 2;
        weight_t sum = Node::get_doublet(c, root->s1, root->s2, 0, 0);
        for (index_t i = 0; i < s; i ++) {
            if (valid.find(i) == valid.end()) continue;
            for (index_t j = 0; j < s; j ++) {
                if (i == j || valid.find(j) == valid.end()) continue;
                graph[0][i][j] = sum - graph[1][i][j];
            }
        }
        c[0] += 1;
        for (index_t i = 0; i < s; i ++) {
            if (valid.find(i) == valid.end()) continue;
            for (index_t j = 0; j < m; j ++) {
                if (valid.find(s + j) == valid.end()) continue;
                sum = Node::get_doublet(c, root->s1, root->s2, j + 1, 0) * root->singlet[j + 1];
                graph[0][s + j][i] = graph[0][i][s + j] = sum - graph[1][i][s + j];
            }
        }
        c[0] += 1;
        for (index_t i = 0; i < m; i ++) {
            if (valid.find(s + i) == valid.end()) continue;
            for (index_t j = 0; j < m; j ++) {
                if (i == j || valid.find(s + j) == valid.end()) continue;
                sum = Node::get_doublet(c, root->s1, root->s2, i + 1, j + 1);
                sum *= root->singlet[i + 1] * root->singlet[j + 1];
                graph[0][s + j][s + i] = graph[0][s + i][s + j] = sum - graph[1][s + i][s + j];
            }
        }
        delete [] c;
    }
}

void Tree::get_quartets(std::unordered_map<quartet_t, weight_t> *quartets) {
    std::vector<Node *> leaves;
    get_leaves(root, &leaves);
    get_depth(root, 0);
    index_t idx[4], size = leaves.size();
    for (idx[0] = 0         ; idx[0] < size; idx[0] ++) 
    for (idx[1] = idx[0] + 1; idx[1] < size; idx[1] ++) 
    for (idx[2] = idx[1] + 1; idx[2] < size; idx[2] ++) 
    for (idx[3] = idx[2] + 1; idx[3] < size; idx[3] ++) {
        Node *lowest = NULL, *nodes[4] = {NULL, NULL, NULL, NULL};
        index_t lowest_count = 0;
        for (index_t i = 0, k = 0; i < 4; i ++) {
            for (index_t j = 0; j < i; j ++) {
                Node *a = leaves[idx[i]], *b = leaves[idx[j]];
                while (a->depth > b->depth) a = a->parent;
                while (b->depth > a->depth) b = b->parent;
                while (a != b) {a = a->parent; b = b->parent;}
                if (lowest == NULL || a->depth > lowest->depth) {
                    lowest = a;
                    lowest_count = 0;
                    nodes[0] = leaves[idx[i]];
                    nodes[1] = leaves[idx[j]];
                }
                if (a == lowest) lowest_count ++;
            }
        }
        if (lowest_count != 1) continue;
        for (index_t i = 0; i < 4; i ++) {
            if (leaves[idx[i]] != nodes[0] && leaves[idx[i]] != nodes[1]) 
                nodes[(nodes[2] == NULL ? 2 : 3)] = leaves[idx[i]];
        }
        index_t indices[4];
        for (index_t i = 0; i < 4; i ++) {
            indices[i] = nodes[i]->index;
        }
        quartet_t quartet = join(indices);
        if (quartets->find(quartet) == quartets->end()) 
            (*quartets)[quartet] = 0;
        (*quartets)[quartet] += 1;
    }
}

void Tree::get_wquartets(std::unordered_map<quartet_t, weight_t> *quartets) {
    std::vector<Node *> leaves;
    get_leaves(root, &leaves);
    get_depth(root, 0);
    index_t idx[4], size = leaves.size();
    for (idx[0] = 0         ; idx[0] < size; idx[0] ++) 
    for (idx[1] = idx[0] + 1; idx[1] < size; idx[1] ++) 
    for (idx[2] = idx[1] + 1; idx[2] < size; idx[2] ++) 
    for (idx[3] = idx[2] + 1; idx[3] < size; idx[3] ++) {
        Node *lowest = NULL, *nodes[4] = {NULL, NULL, NULL, NULL};
        index_t lowest_count = 0;
        for (index_t i = 0, k = 0; i < 4; i ++) {
            for (index_t j = 0; j < i; j ++) {
                Node *a = leaves[idx[i]], *b = leaves[idx[j]];
                while (a->depth > b->depth) a = a->parent;
                while (b->depth > a->depth) b = b->parent;
                while (a != b) {a = a->parent; b = b->parent;}
                if (lowest == NULL || a->depth > lowest->depth) {
                    lowest = a;
                    lowest_count = 0;
                    nodes[0] = leaves[idx[i]];
                    nodes[1] = leaves[idx[j]];
                }
                if (a == lowest) lowest_count ++;
            }
        }
        if (lowest_count != 1) continue;
        for (index_t i = 0; i < 4; i ++) {
            if (leaves[idx[i]] != nodes[0] && leaves[idx[i]] != nodes[1]) 
                nodes[(nodes[2] == NULL ? 2 : 3)] = leaves[idx[i]];
        }
        weight_t l = 0, s = 1;
        Node *a = nodes[0], *b = nodes[1];
        while (a->depth > b->depth) {l += a->length; a = a->parent;}
        while (b->depth > a->depth) {l += b->length; b = b->parent;}
        while (a != b) {l += a->length + b->length; a = a->parent; b = b->parent;}
        Node *c = a;
        a = nodes[2], b = nodes[3];
        while (a->depth > b->depth) {l += a->length; a = a->parent;}
        while (b->depth > a->depth) {l += b->length; b = b->parent;}
        while (a != b) {l += a->length + b->length; a = a->parent; b = b->parent;}
        Node *d = a;
        a = c; b = d;
        while (a->depth > b->depth) {s *= 1 - a->support; a = a->parent;}
        while (b->depth > a->depth) {s *= 1 - b->support; b = b->parent;}
        while (a != b) {s *= (1 - a->support) * (1 - b->support); a = a->parent; b = b->parent;}
        if (a == d) {
            weight_t s0 = 1, s1 = 1;
            a = c; b = nodes[2];
            while (a->depth > b->depth) {s0 *= 1 - a->support; a = a->parent;}
            while (b->depth > a->depth) {b = b->parent;}
            while (a != b) {s0 *= (1 - a->support); a = a->parent; b = b->parent;}
            a = c; b = nodes[3];
            while (a->depth > b->depth) {s1 *= 1 - a->support; a = a->parent;}
            while (b->depth > a->depth) {b = b->parent;}
            while (a != b) {s1 *= (1 - a->support); a = a->parent; b = b->parent;}
            //std::cout << s << ' ' << s0 << ' ' << s1 << std::endl;
            if (s0 > s) s = s0;
            if (s1 > s) s = s1;
        }
        else {
            //std::cout << 0 << std::endl;
        }
        weight_t w = exp(- l) * (1 - s);
        index_t indices[4];
        for (index_t i = 0; i < 4; i ++) {
            indices[i] = nodes[i]->index;
        }
        quartet_t quartet = join(indices);
        // std::cout << quartet << ' ' << w << ' ' << l << ' ' << 1 - s << std::endl;
        if (quartets->find(quartet) == quartets->end()) 
            (*quartets)[quartet] = 0;
        (*quartets)[quartet] += w;
    }
}

void Tree::get_wquartets_(std::unordered_map<quartet_t, weight_t> *quartets) {
    std::vector<Node *> leaves;
    get_leaves(root, &leaves);
    get_depth(root, 0);
    index_t idx[4], size = leaves.size();
    for (idx[0] = 0         ; idx[0] < size; idx[0] ++) 
    for (idx[1] = idx[0] + 1; idx[1] < size; idx[1] ++) 
    for (idx[2] = idx[1] + 1; idx[2] < size; idx[2] ++) 
    for (idx[3] = idx[2] + 1; idx[3] < size; idx[3] ++) {
        Node *lowest = NULL, *nodes[4] = {NULL, NULL, NULL, NULL};
        index_t lowest_count = 0;
        for (index_t i = 0, k = 0; i < 4; i ++) {
            for (index_t j = 0; j < i; j ++) {
                Node *a = leaves[idx[i]], *b = leaves[idx[j]];
                while (a->depth > b->depth) a = a->parent;
                while (b->depth > a->depth) b = b->parent;
                while (a != b) {a = a->parent; b = b->parent;}
                if (lowest == NULL || a->depth > lowest->depth) {
                    lowest = a;
                    lowest_count = 0;
                    nodes[0] = leaves[idx[i]];
                    nodes[1] = leaves[idx[j]];
                }
                if (a == lowest) lowest_count ++;
            }
        }
        if (lowest_count != 1) continue;
        for (index_t i = 0; i < 4; i ++) {
            if (leaves[idx[i]] != nodes[0] && leaves[idx[i]] != nodes[1]) 
                nodes[(nodes[2] == NULL ? 2 : 3)] = leaves[idx[i]];
        }
        weight_t l = 0, s = 1;
        Node *a = nodes[0], *b = nodes[1];
        while (a->depth > b->depth) {l += a->length; a = a->parent;}
        while (b->depth > a->depth) {l += b->length; b = b->parent;}
        while (a != b) {l += a->length + b->length; a = a->parent; b = b->parent;}
        Node *c = a;
        a = nodes[2], b = nodes[3];
        while (a->depth > b->depth) {l += a->length; a = a->parent;}
        while (b->depth > a->depth) {l += b->length; b = b->parent;}
        while (a != b) {l += a->length + b->length; a = a->parent; b = b->parent;}
        Node *d = a;
        a = c; b = d;
        while (a->depth > b->depth) {s *= 1 - a->support; a = a->parent;}
        while (b->depth > a->depth) {s *= 1 - b->support; b = b->parent;}
        while (a != b) {s *= (1 - a->support) * (1 - b->support); a = a->parent; b = b->parent;}
        if (a == d) {
            weight_t s0 = 1, s1 = 1;
            a = c; b = nodes[2];
            while (a->depth > b->depth) {s0 *= 1 - a->support; a = a->parent;}
            while (b->depth > a->depth) {b = b->parent;}
            while (a != b) {s0 *= (1 - a->support); a = a->parent; b = b->parent;}
            a = c; b = nodes[3];
            while (a->depth > b->depth) {s1 *= 1 - a->support; a = a->parent;}
            while (b->depth > a->depth) {b = b->parent;}
            while (a != b) {s1 *= (1 - a->support); a = a->parent; b = b->parent;}
            //std::cout << s << ' ' << s0 << ' ' << s1 << std::endl;
            if (s0 > s) s = s0;
            if (s1 > s) s = s1;
        }
        else {
            //std::cout << 0 << std::endl;
        }
        weight_t w = 1 - s;
        index_t indices[4];
        for (index_t i = 0; i < 4; i ++) {
            indices[i] = nodes[i]->index;
        }
        quartet_t quartet = join(indices);
        // std::cout << quartet << ' ' << w << ' ' << l << ' ' << 1 - s << std::endl;
        if (quartets->find(quartet) == quartets->end()) 
            (*quartets)[quartet] = 0;
        (*quartets)[quartet] += w;
    }
}

std::string Tree::to_string(std::unordered_map<quartet_t, weight_t> &quartets) {
    std::string s = "";
    for (auto elem : quartets) {
        index_t *indices = split(elem.first);
        s += dict->index2label(indices[0]) + "," + dict->index2label(indices[1]) + "|";
        s += dict->index2label(indices[2]) + "," + dict->index2label(indices[3]) + ":";
        s += std::to_string((double)elem.second) + ";\n";
        delete [] indices;
    }
    return s;
}

void Tree::get_leaves(Node *root, std::vector<Node *> *leaves) {
    if (root->children.size() == 0) 
        leaves->push_back(root);
    for (Node *child : root->children) 
        get_leaves(child, leaves);
}

void Tree::get_leaf_set(Node *root, std::unordered_set<Node *> *leaf_set) {
    if (root->children.size() == 0) 
        leaf_set->insert(root);
    for (Node *child : root->children)
        get_leaf_set(child, leaf_set);
}

void Tree::get_depth(Node *root, index_t depth) {
    root->depth = depth;
    for (Node *child : root->children) 
        get_depth(child, depth + 1);
}

Node* Tree::find_node_for_split(std::unordered_set<index_t> &clade) {
    std::vector<Node*> nodes, postorder_nodes;
    std::queue<Node*> queue;
    Node *node;
    index_t ni, ci, n_leaves = 0, n_nodes = 0;
    std::unordered_map<index_t, index_t> ugh;

    // Prepare for postorder traversal and count
    nodes.push_back(root);
    while (nodes.size() > 0) {
        node = nodes.back();
        nodes.pop_back();
        for (Node* child : node->children)
            nodes.push_back(child);
        postorder_nodes.push_back(node);

        if (node->is_leaf()) n_leaves++;

	ugh.insert({node->index, n_nodes});
        n_nodes++;
    }

    // Run postorder traversal and count
    std::vector<index_t> n_lv_below(n_nodes, 0);
    std::vector<index_t> n_og_below(n_nodes, 0);

    while (postorder_nodes.size() > 0) {
        node = postorder_nodes.back();
        postorder_nodes.pop_back();

        ni = ugh[node->index];

        // Update counts
        if (node->is_leaf()) {
            if (clade.find(node->index) != clade.end())
                n_og_below[ni]++;
            n_lv_below[ni]++;
        } else {
            for (Node* child : node->children) {
                ci = ugh[child->index];
                n_og_below[ni] += n_og_below[ci];
                n_lv_below[ni] += n_lv_below[ci];
            }
        }

        // Check if found clade
        if (n_og_below[ni] == clade.size() && n_lv_below[ni] == clade.size()) {
            return node;
        }
    }

    // Run preorder traversal
    queue.push(root);
    while (queue.size() > 0) {
        node = queue.front();
        queue.pop();
        for (Node* child: node->children)
            queue.push(child);

        // Check if clade is split across root
        ni = ugh[node->index];
        if (n_og_below[ni] == 0) {
            if ((n_leaves - n_lv_below[ni]) == clade.size()) {
                return node;
            }
        }
    }

    return NULL;
}

Node* Tree::find_node(index_t index) {
    std::queue<Node*> queue;
    Node *node;

    // Run preorder traversal
    queue.push(root);
    while (queue.size() > 0) {
        node = queue.front();
        queue.pop();

        if (node->index == index)
            return node;

        for (Node* child: node->children)
            queue.push(child);
    }

    return NULL;
}

void Tree::reroot_on_edge_above_node(Node *node) {
    if (node == NULL) return;

    std::stack<Node*> stack;
    Node *old_root, *next_root, *sibl, *keep;

    keep = node;

    while (node != NULL) {
        stack.push(node);
        node = node->get_parent();
    }
    stack.pop();  // Remove root from stack

    // Start relocating root
    // Basically swapping nodes on stack across root
    node = stack.top();
    stack.pop();
    while (node != keep) {
        // Move branch length over root
        sibl = node->get_sibling();
        sibl->length += node->length;

        // Get next root
        next_root = stack.top();

        // Remove node from root
        old_root = node->get_parent();
        old_root->remove_child(node);

        // Remove next root's sibling 
        sibl = next_root->get_sibling();
        node->remove_child(sibl);

        // Attach next root's sibling to old root
        old_root->add_child(sibl);
        old_root->f[0] = next_root->f[0]; // need to do for all data
        old_root->f[1] = next_root->f[1];
        old_root->f[2] = next_root->f[2];
        old_root->support = next_root->support;
        if (next_root == pcs_node) pcs_node = old_root;

        // Make node the new root so next root is below it
        node->add_child(old_root);
        node->f[0] = 0; // need to do for all data
        node->f[1] = 0;
        node->f[2] = 0;
        node->support = 0.0;
        node->length = 0.0;
        this->root = node;

        node = next_root;
        stack.pop();
    }
}
