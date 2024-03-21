#include "tree.hpp"

void Tree::build_wstates(Node *root) {
    root->size = 4;
    root->ssinglet = init(root->size + 1);
    root->ssinglet_ = init(root->size + 1);
    for (index_t k = 0; k < 2; k ++) {
        root->mdoublet[k] = init(root->size + 1);
        root->mdoublet_[k] = init(root->size + 1);
        root->sdoublet[k] = new weight_t*[root->size + 1];
        for (index_t i = 0; i <= root->size; i ++) 
            root->sdoublet[k][i] = init(root->size + 1);
        root->sdoublet_[k] = new weight_t*[root->size + 1];
        for (index_t i = 0; i <= root->size; i ++) 
            root->sdoublet_[k][i] = init(root->size + 1);
    }
    for (Node *child : root->children) 
        build_wstates(child);
}

void Tree::clear_wstates_(Node *root) {
    delete [] root->ssinglet;
    delete [] root->ssinglet_;
    for (index_t k = 0; k < 2; k ++) {
        delete [] root->mdoublet[k];
        delete [] root->mdoublet_[k];
        for (index_t i = 0; i <= root->size; i ++) 
            delete [] root->sdoublet[k][i];
        delete [] root->sdoublet[k];
        for (index_t i = 0; i <= root->size; i ++) 
            delete [] root->sdoublet_[k][i];
        delete [] root->sdoublet_[k];
    }
    for (Node *child : root->children) 
        clear_wstates_(child);
}

void Tree::build_ssinglet(Node *root, std::unordered_map<index_t, index_t> quad) {
    if (root->children.size() == 0) {
        for (index_t i = 0; i <= root->size; i ++)
            root->ssinglet[i] = 0;
        root->ssinglet[quad[root->index]] ++;
    }
    else {
        for (Node *child : root->children) 
            build_ssinglet(child, quad);
        for (index_t i = 0; i <= root->size; i ++) {
            root->ssinglet[i] = 
                root->children[0]->ssinglet[i] * root->children[0]->length_ + 
                root->children[1]->ssinglet[i] * root->children[1]->length_;
        }
    }
}

weight_t Tree::get_qcount(std::unordered_map<index_t, index_t> quad) {
    index_t c1 = 0;
    index_t c2 = 0;
    index_t c3 = 0;
    index_t c4 = 0;
    for ( const auto &myPair : quad ) {
        Node* nodeptr = index2node[myPair.first];

        if (nodeptr != NULL) {
            if (myPair.second == 1) c1++;
            else if (myPair.second == 2) c2++;
            else if (myPair.second == 3) c3++;
            else c4++;
        }
    }
    return c1 * c2 * c3 * c4;
}

weight_t Tree::get_qfreq(std::unordered_map<index_t, index_t> quad) {
    build_wstates(root);
    build_ssinglet(root, quad);
    build_ssinglet_(root);
    build_sdoublet(root);
    build_sdoublet_(root);
    weight_t s = freq_(root);
    clear_wstates_(root);
    return s;
}

weight_t Tree::freq_(Node *root) {
    weight_t f = 0;
    if (root->children.size() != 0) {
        for (Node *child : root->children) 
            f += freq_(child);
        Node *left = root->children[0], *right = root->children[1];
        weight_t s[2] = {0, 0};
        for (index_t k = 0; k < 2; k ++) {
            s[k] += left->ssinglet[1] * left->length_ * right->ssinglet[2] * right->length_ * root->sdoublet_[k][3][4];
            s[k] += left->ssinglet[2] * left->length_ * right->ssinglet[1] * right->length_ * root->sdoublet_[k][3][4];
            s[k] += left->sdoublet[k][3][4] * left->support_[k] * right->ssinglet[1] * right->length_ * root->ssinglet_[2];
            s[k] += left->sdoublet[k][3][4] * left->support_[k] * right->ssinglet[2] * right->length_ * root->ssinglet_[1];
            s[k] += left->ssinglet[1] * left->length_ * right->sdoublet[k][3][4] * right->support_[k] * root->ssinglet_[2];
            s[k] += left->ssinglet[2] * left->length_ * right->sdoublet[k][3][4] * right->support_[k] * root->ssinglet_[1];
            s[k] += left->sdoublet[k][1][2] * left->support_[k] * right->ssinglet[3] * right->length_ * root->ssinglet_[4];
            s[k] += left->sdoublet[k][1][2] * left->support_[k] * right->ssinglet[4] * right->length_ * root->ssinglet_[3];
            s[k] += left->ssinglet[3] * left->length_ * right->sdoublet[k][1][2] * right->support_[k] * root->ssinglet_[4];
            s[k] += left->ssinglet[4] * left->length_ * right->sdoublet[k][1][2] * right->support_[k] * root->ssinglet_[3];
            s[k] += left->ssinglet[3] * left->length_ * right->ssinglet[4] * right->length_ * root->sdoublet_[k][1][2];
            s[k] += left->ssinglet[4] * left->length_ * right->ssinglet[3] * right->length_ * root->sdoublet_[k][1][2];
        }
        f += s[1] - s[0];
    }
    return f;
}

void Tree::build_wstates_s(Node *root) {
    root->size = 0;
    root->ssinglet = init(root->size + 1);
    root->ssinglet_ = init(root->size + 1);
    for (index_t k = 0; k < 2; k ++) {
        root->mdoublet[k] = init(root->size + 1);
        root->mdoublet_[k] = init(root->size + 1);
        root->sdoublet[k] = new weight_t*[root->size + 1];
        for (index_t i = 0; i <= root->size; i ++) 
            root->sdoublet[k][i] = init(root->size + 1);
        root->sdoublet_[k] = new weight_t*[root->size + 1];
        for (index_t i = 0; i <= root->size; i ++) 
            root->sdoublet_[k][i] = init(root->size + 1);
    }
    for (Node *child : root->children) 
        build_wstates_s(child);
}

void Tree::build_ssinglet_s(Node *root) {
    if (root->children.size() == 0) {
        for (index_t i = 0; i <= root->size; i ++)
            root->ssinglet[i] = 0;
        root->ssinglet[0] ++;
    }
    else {
        for (Node *child : root->children) 
            build_ssinglet_s(child);
        for (index_t i = 0; i <= root->size; i ++) {
            root->ssinglet[i] = 
                root->children[0]->ssinglet[i] * root->children[0]->length_ + 
                root->children[1]->ssinglet[i] * root->children[1]->length_;
        }
    }
}

weight_t Tree::total_weight() {
    build_wstates_s(root);
    build_ssinglet_s(root);
    build_ssinglet_(root);
    build_sdoublet(root);
    build_sdoublet_(root);
    total_quartet_weight = freq_s(root) / 2;
    // std::cout << "test qs: " << total_quartet_weight << ' ' << total_weight_bf() << std::endl;
    clear_wstates_(root);
    return total_quartet_weight;
}

weight_t Tree::freq_s(Node *root) {
    weight_t f = 0;
    if (root->children.size() != 0) {
        for (Node *child : root->children) 
            f += freq_s(child);
        Node *left = root->children[0], *right = root->children[1];
        weight_t s[2] = {0, 0};
        for (index_t k = 0; k < 2; k ++) {
            s[k] += left->ssinglet[0] * left->length_ * right->ssinglet[0] * right->length_ * root->sdoublet_[k][0][0];
            s[k] += left->sdoublet[k][0][0] * left->support_[k] * right->ssinglet[0] * right->length_ * root->ssinglet_[0];
            s[k] += left->ssinglet[0] * left->length_ * right->sdoublet[k][0][0] * right->support_[k] * root->ssinglet_[0];
        }
        f += s[1] - s[0];
    }
    return f;
}

weight_t Tree::total_weight_bf() {
    std::unordered_map<quartet_t, weight_t> q;
    get_wquartets(&q);
    weight_t s = 0;
    for (auto elem : q) {
        s += elem.second;
    }
    return s;
}
