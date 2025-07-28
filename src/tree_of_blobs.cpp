#include "tree.hpp"

// O(n*cn*klogn), O(kn^3logn) if c is O(n) 

SpeciesTree::SpeciesTree(Tree *input, Dict *dict, weight_t alpha, weight_t beta) {
    std::cout << "Constructing tree of blobs from annotated species tree" << std::endl;
    this->dict = dict;
    std::vector<Node *> internal;
    std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> bips;
    input->get_bipartitions(&internal, &bips);
    std::cout << bips.size() << " branches to test" << std::endl;
    std::unordered_set<Node *> false_positive;
    for (index_t i = 0; i < internal.size(); i ++) {
        if (internal[i]->min_pvalue < alpha || internal[i]->max_pvalue > beta) 
            false_positive.insert(internal[i]);
    }
    root = build_refinement(input->root, false_positive);
}

SpeciesTree::SpeciesTree(std::vector<Tree *> &input, Dict *dict, SpeciesTree* display, weight_t alpha, weight_t beta, unsigned long int iter_limit_blob) {
    std::cout << "Constructing tree of blobs" << std::endl;
    RINS.parseEvalQ("library(MSCquartets, lib.loc=\"/fs/cbcb-lab/ekmolloy/yhhan/tree-of-blobs/software/Rlibs\")");
    for (Tree *t : input) t->LCA_preprocessing();
    display->refine();
    std::string mode = "n";
    display->annotate(input, mode);
    std::vector<Node *> internal;
    std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> bips;
    display->get_bipartitions(&internal, &bips);
    std::cout << bips.size() << " branches to test" << std::endl;
    std::unordered_set<Node *> false_positive;
    size_t iter_limit = iter_limit_blob;
    for (index_t i = 0; i < internal.size(); i ++) {
        std::cout << "Testing branch " << i << ", ";
        if ((alpha >= 0 || beta <= 1) && internal[i]->isfake) {
            false_positive.insert(internal[i]);
            std::cout << "fake ***" << std::endl;
            continue;
        }
        weight_t min, max;
        if (iter_limit != 0) {
            min = search(input, bips[i].first, bips[i].second, iter_limit, internal[i]->min_f);
            // max = search_star(input, bips[i].first, bips[i].second, iter_limit, internal[i]->max_f);
        }
	else {
            min = search(input, bips[i].first, bips[i].second, internal[i]->min_f);
            // max = search_star(input, bips[i].first, bips[i].second, internal[i]->max_f);
	}
        max = pvalue_star(internal[i]->f);
        internal[i]->min_pvalue = min;
        internal[i]->max_pvalue = max;
        std::cout << "QTT: " << min << " ";
        std::cout << "[" << internal[i]->min_f[0] << "/" << internal[i]->min_f[1] << "/" << internal[i]->min_f[2] << "] ";
 	std::cout << "QST: " << max << " ";
        std::cout << "[" << internal[i]->f[0] << "/" << internal[i]->f[1] << "/" << internal[i]->f[2] << "]";
        if (min < alpha || max > beta) {
            false_positive.insert(internal[i]);
            std::cout << " ***";
        }
        std::cout << std::endl;
    }
    if (display->root->children.size() == 2) 
        false_positive.insert(display->root->children[1]);
    this->dict = display->dict;
    root = build_refinement(display->root, false_positive);
}

Node *SpeciesTree::build_refinement(Node *root, std::unordered_set<Node *> false_positive) {
    if (root->children.size() == 0) {
        Node *new_root = new Node(root->index, false);
        index2node[new_root->index] = new_root;
        return new_root;
    }
    else {
        Node *new_root = new Node(pseudonym(), false_positive.find(root) != false_positive.end());
        for (Node *child : root->children) {
            Node *new_child = build_refinement(child, false_positive);
            if (! new_child->isfake) {
                new_root->children.push_back(new_child);
            }
            else {
                for (Node *grand_child : new_child->children) 
                    new_root->children.push_back(grand_child);
                new_child->children.clear();
                delete new_child;
            }
        }
        for (Node *new_child : new_root->children) 
            new_child->parent = new_root;
        return new_root;
    }
}

void Tree::get_bipartition(Node *root, std::vector<Node *> *A, std::vector<Node *> *B) {
    std::vector<Node *> leaves;
    get_leaves(this->root, &leaves);
    std::unordered_set<Node *> leaf_set;
    get_leaf_set(root, &leaf_set);
    for (Node *leaf : leaves) {
        if (leaf_set.find(leaf) == leaf_set.end()) 
            A->push_back(leaf);
        else 
            B->push_back(leaf);
    }
}

void Tree::get_bipartitions(Node *root, std::vector<Node *> *internal, std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> *bips) {
    if (root->children.size() == 0) return ;
    if (root->parent != NULL && (root->parent->parent != NULL || root->parent->children.size() > 2 || root == root->parent->children[0])) {
        std::vector<Node *> A, B;
        get_bipartition(root, &A, &B);
        bips->push_back(std::make_pair(A, B));
        internal->push_back(root); // std::cout << "?" << root->support << std::endl;
    }
    for (Node *child : root->children) {
        get_bipartitions(child, internal, bips);
    }
}

void Tree::get_bipartitions(std::vector<Node *> *internal, std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> *bips) {
    get_bipartitions(root, internal, bips);
}

std::string Tree::display_bipartition(std::vector<Node *> &A, std::vector<Node *> &B) {
    std::string s = "";
    for (Node *a : A) s += std::to_string(a->index) + " ";
    s += "| ";
    for (Node *b : B) s += std::to_string(b->index) + " ";
    return s;
}

weight_t SpeciesTree::search(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B, weight_t *min_f) {
    index_t i[4];
    weight_t min = -1;
    size_t count = 0;
    for (i[0] = 0; i[0] < A.size(); i[0] ++) {
        for (i[1] = i[0] + 1; i[1] < A.size(); i[1] ++) {
            for (i[2] = 0; i[2] < B.size(); i[2] ++) {
                for (i[3] = i[2] + 1; i[3] < B.size(); i[3] ++) {
                    index_t temp[4];
                    temp[0] = A[i[0]]->index;
                    temp[1] = A[i[1]]->index;
                    temp[2] = B[i[2]]->index;
                    temp[3] = B[i[3]]->index;
                    count ++;
                    weight_t f[3];
                    weight_t score = get_pvalue(input, temp, f);
                    if (min < 0 || score < min) {
                        min = score;
                        min_f[0] = f[0]; min_f[1] = f[1]; min_f[2] = f[2];
                    }
                }
            }
        }
    }
    //std::cout << "naive iter: " << count << std::endl;
    return min;
}

weight_t SpeciesTree::search_star(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B, weight_t *min_f) {
    index_t i[4];
    weight_t min = -1;
    size_t count = 0;
    for (i[0] = 0; i[0] < A.size(); i[0] ++) {
        for (i[1] = i[0] + 1; i[1] < A.size(); i[1] ++) {
            for (i[2] = 0; i[2] < B.size(); i[2] ++) {
                for (i[3] = i[2] + 1; i[3] < B.size(); i[3] ++) {
                    index_t temp[4];
                    temp[0] = A[i[0]]->index;
                    temp[1] = A[i[1]]->index;
                    temp[2] = B[i[2]]->index;
                    temp[3] = B[i[3]]->index;
                    count ++;
                    weight_t f[3];
                    weight_t score = get_pvalue_star(input, temp, f);
                    if (min < 0 || score > min) {
                        min = score;
                        min_f[0] = f[0]; min_f[1] = f[1]; min_f[2] = f[2];
                    }
                }
            }
        }
    }
    //std::cout << "naive iter: " << count << std::endl;
    return min;
}

weight_t SpeciesTree::search(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B, size_t iter_limit, weight_t *min_f) {
    index_t indices[4];
    weight_t min = -1;
    size_t count = 0;
    while (count < iter_limit) {
        indices[0] = A[rand() % A.size()]->index;
        do {indices[1] = A[rand() % A.size()]->index;} while (indices[0] == indices[1]);
        indices[2] = B[rand() % B.size()]->index;
        do {indices[3] = B[rand() % B.size()]->index;} while (indices[2] == indices[3]);
        count += neighbor_search(input, A, B, indices, &min, min_f);
        // std::cout << i << ' ' << min << std::endl;
    }
    //std::cout << "heuristic iter: " << count << std::endl;
    return min;
}

weight_t SpeciesTree::neighbor_search(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B, index_t *current, weight_t *min, weight_t *min_f) {
    weight_t current_f[3];
    weight_t current_score = get_pvalue(input, current, current_f);
    size_t k = 1;
    while (true) {
        index_t temp[4], next[4];
        weight_t temp_score, next_score = -1, temp_f[3], next_f[3];
        for (index_t j = 0; j < 4; j ++) 
            temp[j] = current[j];
        for (index_t i = 0; i < A.size(); i ++) {
            index_t new_index = A[i]->index;
            if (new_index == current[0] || new_index == current[1]) continue;
            temp[0] = new_index;
            temp_score = get_pvalue(input, temp, temp_f); k ++;
            if (next_score < 0 || temp_score < next_score) {
                next_score = temp_score;
                for (index_t j = 0; j < 4; j ++) 
                    next[j] = temp[j];
                next_f[0] = temp_f[0]; next_f[1] = temp_f[1]; next_f[2] = temp_f[2];
            }
            temp[0] = current[0];
            temp[1] = new_index;
            temp_score = get_pvalue(input, temp, temp_f); k ++;
            if (next_score < 0 || temp_score < next_score) {
                next_score = temp_score;
                for (index_t j = 0; j < 4; j ++) 
                    next[j] = temp[j];
                next_f[0] = temp_f[0]; next_f[1] = temp_f[1]; next_f[2] = temp_f[2];
            }
            temp[1] = current[1];
        }
        for (index_t i = 0; i < B.size(); i ++) {
            index_t new_index = B[i]->index;
            if (new_index == current[2] || new_index == current[3]) continue;
            temp[2] = new_index;
            temp_score = get_pvalue(input, temp, temp_f); k ++;
            if (next_score < 0 || temp_score < next_score) {
                next_score = temp_score;
                for (index_t j = 0; j < 4; j ++) 
                    next[j] = temp[j];
                next_f[0] = temp_f[0]; next_f[1] = temp_f[1]; next_f[2] = temp_f[2];
            }
            temp[2] = current[2];
            temp[3] = new_index;
            temp_score = get_pvalue(input, temp, temp_f); k ++;
            if (next_score < 0 || temp_score < next_score) {
                next_score = temp_score;
                for (index_t j = 0; j < 4; j ++) 
                    next[j] = temp[j];
                next_f[0] = temp_f[0]; next_f[1] = temp_f[1]; next_f[2] = temp_f[2];
            }
            temp[3] = current[3];
        }
        if (next_score >= current_score) break;
        current_score = next_score;
        for (index_t j = 0; j < 4; j ++) 
            current[j] = next[j];
        current_f[0] = next_f[0]; current_f[1] = next_f[1]; current_f[2] = next_f[2];
    }
    if (*min < 0 || *min > current_score) {
        *min = current_score;
        min_f[0] = current_f[0]; min_f[1] = current_f[1]; min_f[2] = current_f[2];
    }
    return k;
}

weight_t SpeciesTree::search_star(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B, size_t iter_limit, weight_t *min_f) {
    index_t indices[4];
    weight_t min = -1;
    size_t count = 0;
    while (count < iter_limit) {
        indices[0] = A[rand() % A.size()]->index;
        do {indices[1] = A[rand() % A.size()]->index;} while (indices[0] == indices[1]);
        indices[2] = B[rand() % B.size()]->index;
        do {indices[3] = B[rand() % B.size()]->index;} while (indices[2] == indices[3]);
        count += neighbor_search_star(input, A, B, indices, &min, min_f);
        // std::cout << i << ' ' << min << std::endl;
    }
    //std::cout << "heuristic iter: " << count << std::endl;
    return min;
}

weight_t SpeciesTree::neighbor_search_star(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B, index_t *current, weight_t *min, weight_t *min_f) {
    weight_t current_f[3];
    weight_t current_score = get_pvalue_star(input, current, current_f);
    size_t k = 1;
    while (true) {
        index_t temp[4], next[4];
        weight_t temp_score, next_score = -1, temp_f[3], next_f[3];
        for (index_t j = 0; j < 4; j ++) 
            temp[j] = current[j];
        for (index_t i = 0; i < A.size(); i ++) {
            index_t new_index = A[i]->index;
            if (new_index == current[0] || new_index == current[1]) continue;
            temp[0] = new_index;
            temp_score = get_pvalue_star(input, temp, temp_f); k ++;
            if (next_score < 0 || temp_score > next_score) {
                next_score = temp_score;
                for (index_t j = 0; j < 4; j ++) 
                    next[j] = temp[j];
                next_f[0] = temp_f[0]; next_f[1] = temp_f[1]; next_f[2] = temp_f[2];
            }
            temp[0] = current[0];
            temp[1] = new_index;
            temp_score = get_pvalue_star(input, temp, temp_f); k ++;
            if (next_score < 0 || temp_score > next_score) {
                next_score = temp_score;
                for (index_t j = 0; j < 4; j ++) 
                    next[j] = temp[j];
                next_f[0] = temp_f[0]; next_f[1] = temp_f[1]; next_f[2] = temp_f[2];
            }
            temp[1] = current[1];
        }
        for (index_t i = 0; i < B.size(); i ++) {
            index_t new_index = B[i]->index;
            if (new_index == current[2] || new_index == current[3]) continue;
            temp[2] = new_index;
            temp_score = get_pvalue_star(input, temp, temp_f); k ++;
            if (next_score < 0 || temp_score > next_score) {
                next_score = temp_score;
                for (index_t j = 0; j < 4; j ++) 
                    next[j] = temp[j];
                next_f[0] = temp_f[0]; next_f[1] = temp_f[1]; next_f[2] = temp_f[2];
            }
            temp[2] = current[2];
            temp[3] = new_index;
            temp_score = get_pvalue_star(input, temp, temp_f); k ++;
            if (next_score < 0 || temp_score > next_score) {
                next_score = temp_score;
                for (index_t j = 0; j < 4; j ++) 
                    next[j] = temp[j];
                next_f[0] = temp_f[0]; next_f[1] = temp_f[1]; next_f[2] = temp_f[2];
            }
            temp[3] = current[3];
        }
        if (next_score <= current_score) break;
        current_score = next_score;
        for (index_t j = 0; j < 4; j ++) 
            current[j] = next[j];
        current_f[0] = next_f[0]; current_f[1] = next_f[1]; current_f[2] = next_f[2];
    }
    if (*min < 0 || *min < current_score) {
        *min = current_score;
        min_f[0] = current_f[0]; min_f[1] = current_f[1]; min_f[2] = current_f[2];
    }
    return k;
}

weight_t SpeciesTree::get_pvalue(std::vector<Tree *> &input, index_t *indices, weight_t *f) {
    index_t temp[4];
    for (index_t i = 0; i < 4; i ++) 
        temp[i] = indices[i];
    std::sort(temp, temp + 4);
    quartet_t q = join(temp);
    if (pvalues.find(q) == pvalues.end()) {
        weight_t qCF[3] = {0, 0, 0};
        for (Tree *t : input) {
            index_t topology = t->get_quartet(temp);
            if (topology >= 0) qCF[topology] += 1;
        }
        f[0] = qCF[0]; f[1] = qCF[1]; f[2] = qCF[2];
        /*
        std::sort(qCF, qCF + 3);
        std::vector<weight_t> all = pvalue_all(indices);
        std::sort(all.begin(), all.begin() + 3);
        bool verification = qCF[0] == all[0] && qCF[1] == all[1] && qCF[2] == all[2];
        if (! verification) {
            std::cout << qCF[0] << ' ' << qCF[1] << ' ' << qCF[2] << std::endl;
            std::cout << all[0] << ' ' << all[1] << ' ' << all[2] << std::endl;
        }
        assert(verification);
        */
        // std::cout << qCF[0] << ' ' << qCF[1] << ' ' << qCF[2] << std::endl;
	pvalues[q] = std::make_pair(pvalue(qCF), f);
        // pvalues[q] = pvalue(indices);

    }
    auto elem = pvalues[q];
    f[0] = elem.second[0]; f[1] = elem.second[1]; f[2] = elem.second[2];
    return elem.first;
}

weight_t SpeciesTree::get_pvalue_star(std::vector<Tree *> &input, index_t *indices, weight_t *f) {
    index_t temp[4];
    for (index_t i = 0; i < 4; i ++) 
        temp[i] = indices[i];
    std::sort(temp, temp + 4);
    quartet_t q = join(temp);
    if (pvalues_star.find(q) == pvalues_star.end()) {
        weight_t qCF[3] = {0, 0, 0};
        for (Tree *t : input) {
            index_t topology = t->get_quartet(temp);
            if (topology >= 0) qCF[topology] += 1;
        }
        f[0] = qCF[0]; f[1] = qCF[1]; f[2] = qCF[2];
        /*
        std::sort(qCF, qCF + 3);
        std::vector<weight_t> all = pvalue_all(indices);
        std::sort(all.begin(), all.begin() + 3);
        bool verification = qCF[0] == all[0] && qCF[1] == all[1] && qCF[2] == all[2];
        if (! verification) {
            std::cout << qCF[0] << ' ' << qCF[1] << ' ' << qCF[2] << std::endl;
            std::cout << all[0] << ' ' << all[1] << ' ' << all[2] << std::endl;
        }
        assert(verification);
        */
        // std::cout << qCF[0] << ' ' << qCF[1] << ' ' << qCF[2] << std::endl;
	pvalues_star[q] = std::make_pair(pvalue_star(qCF), f);
        // pvalues[q] = pvalue(indices);

    }
    auto elem = pvalues_star[q];
    f[0] = elem.second[0]; f[1] = elem.second[1]; f[2] = elem.second[2];
    return elem.first;
}

index_t Tree::get_quartet(index_t *indices) {
    Node *leaves[4];
    for (index_t i = 0; i < 4; i ++) {
        leaves[i] = index2node[indices[i]];
    }
    Node *lowest = NULL;
    index_t nodes[4] = {-1, -1, -1, -1};
    index_t lowest_count = 0;
    for (index_t i = 0; i < 4; i ++) {
        for (index_t j = i + 1; j < 4; j ++) {
            // Node *a = LCA_naive(leaves[i], leaves[j]);
            Node *a = LCA_fast(leaves[i], leaves[j]);
            if (lowest == NULL || a->depth > lowest->depth) {
                lowest = a;
                lowest_count = 0;
                nodes[0] = i;
                nodes[1] = j;
            }
            if (a == lowest) lowest_count ++;
        }
    }
    if (lowest_count != 1) return -1;
    if (nodes[0] == 0) 
        return nodes[1] - 1;
    else 
        return 5 - nodes[0] - nodes[1];
}

void Tree::LCA_preprocessing() {
    std::vector<Node *> stack;
    stack.push_back(root);
    LCA_depth_first_search(root, stack);
}

void Tree::LCA_depth_first_search(Node *root, std::vector<Node *> &stack) {
    root->depth = stack.size() - 1;
    for (index_t i = 1; i < stack.size(); i *= 2) 
        root->ancestors.push_back(stack[stack.size() - i - 1]);
    if (root->children.size() != 0) {
        for (Node *child : root->children) {
            stack.push_back(child);
            LCA_depth_first_search(child, stack);
            stack.pop_back();
        }
    }
}

Node *Tree::LCA_fast(Node *x, Node *y) {
    if (x->depth > y->depth) {
        Node *temp = x;
        x = y;
        y = temp;
    }
    for (index_t i = y->ancestors.size(); i >= 0; i --) {
        if (i >= y->ancestors.size()) continue;
        if (y->ancestors[i]->depth >= x->depth) 
            y = y->ancestors[i];
    }
    if (x == y) return y;
    for (index_t i = y->ancestors.size(); i >= 0; i --) {
        if (i >= y->ancestors.size()) continue;
        if (y->ancestors[i] != x->ancestors[i]) {
            y = y->ancestors[i];
            x = x->ancestors[i];
        }
    }
    return y->ancestors[0];
}

Node *Tree::LCA_naive(Node *a, Node *b) {
    while (a->depth > b->depth) a = a->parent;
    while (b->depth > a->depth) b = b->parent;
    while (a != b) {a = a->parent; b = b->parent;}
    return a;
}
