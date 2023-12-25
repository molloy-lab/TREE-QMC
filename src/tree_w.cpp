#include "tree.hpp"

extern bool DEBUG_MODE;

void Tree::clear_wstates(Node *root) {
    delete [] root->ssinglet;
    delete [] root->ssinglet_;
    for (index_t k = 0; k < 2; k ++) {
        delete [] root->pdoublet[k];
        delete [] root->ptriplet[k];
        delete [] root->mdoublet[k];
        delete [] root->mdoublet_[k];
        for (index_t i = 0; i <= root->size; i ++) 
            delete [] root->sdoublet[k][i];
        delete [] root->sdoublet[k];
        for (index_t i = 0; i <= root->size; i ++) 
            delete [] root->sdoublet_[k][i];
        delete [] root->sdoublet_[k];
        for (index_t i = 0; i <= root->size; i ++) 
            delete [] root->striplet[k][i];
        delete [] root->striplet[k];
    }
    for (Node *child : root->children) 
        clear_wstates(child);
}

void Tree::build_wstates(Node *root, Taxa &subset) {
    root->plength = std::nan("");
    root->tdoublet[0] = root->tdoublet[1] = std::nan("");
    root->tdoublet_[0] = root->tdoublet_[1] = std::nan("");
    root->size = subset.artificial_taxa();
    root->ssinglet = init(root->size + 1);
    root->ssinglet_ = init(root->size + 1);
    for (index_t k = 0; k < 2; k ++) {
        root->pdoublet[k] = init(root->size + 1);
        root->ptriplet[k] = init(root->size + 1);
        root->mdoublet[k] = init(root->size + 1);
        root->mdoublet_[k] = init(root->size + 1);
        root->sdoublet[k] = new weight_t*[root->size + 1];
        for (index_t i = 0; i <= root->size; i ++) 
            root->sdoublet[k][i] = init(root->size + 1);
        root->sdoublet_[k] = new weight_t*[root->size + 1];
        for (index_t i = 0; i <= root->size; i ++) 
            root->sdoublet_[k][i] = init(root->size + 1);
        root->striplet[k] = new weight_t*[root->size + 1];
        for (index_t i = 0; i <= root->size; i ++) 
            root->striplet[k][i] = init(root->size + 1);
    }
    for (Node *child : root->children) 
        build_wstates(child, subset);
}

void Tree::build_ssinglet(Node *root, Taxa &subset) {
    if (root->children.size() == 0) {
        for (index_t i = 0; i <= root->size; i ++)
            root->ssinglet[i] = 0;
        index_t index = subset.root_key(root->index);
        root->ssinglet[index] = subset.root_weight(root->index);
    }
    else {
        for (Node *child : root->children) 
            build_ssinglet(child, subset);
        for (index_t i = 0; i <= root->size; i ++) {
            root->ssinglet[i] = 
                root->children[0]->ssinglet[i] * root->children[0]->length_ + 
                root->children[1]->ssinglet[i] * root->children[1]->length_;
        }
    }
}

void Tree::build_ssinglet_(Node *root) {
    if (root->parent == NULL) {
        for (index_t i = 0; i <= root->size; i ++)
            root->ssinglet_[i] = 0;
    }
    if (root->children.size() == 0) return ;
    for (index_t c = 0; c < 2; c ++) {
        for (index_t i = 0; i <= root->size; i ++) {
            root->children[c]->ssinglet_[i] = 
                root->ssinglet_[i] * root->children[c]->length_ + 
                root->children[1 - c]->ssinglet[i] * root->children[1 - c]->length_ * root->children[c]->length_;
        }
        build_ssinglet_(root->children[c]);
    }
}

template <typename function1, typename function2>
weight_t Tree::squartet(function1 f1, function2 f2, index_t size, index_t x, index_t y) {
    weight_t s1 = 0, s2 = 0, s12 = 0;
    for (index_t i = 0; i <= root->size; i ++) {
        if (i != 0 && (i == x || i == y)) continue;
        weight_t w1 = f1(i), w2 = f2(i);
        s1 += w1; s2 += w2;
        if (i != 0) s12 += w1 * w2;
    }
    return s1 * s2 - s12;
}

void Tree::build_sdoublet(Node *root) {
    if (root->children.size() == 0) {
        for (index_t i = 0; i <= root->size; i ++) {
            for (index_t j = 0; j <= root->size; j ++) 
                root->sdoublet[0][i][j] = root->sdoublet[1][i][j] = 0;
            root->mdoublet[0][i] = root->mdoublet[1][i] = 0;
        }
        root->tdoublet[0] = root->tdoublet[1] = 0;
    }
    else {
        for (Node *child : root->children) 
            build_sdoublet(child);
        root->tdoublet[0] = root->tdoublet[1] = 0;
        for (index_t i = 0; i <= root->size; i ++) {
            root->mdoublet[0][i] = root->mdoublet[1][i] = 0;
            for (index_t j = 0; j <= root->size; j ++) {
                if (i > 0 && i == j) {
                    root->sdoublet[0][i][j] = 0;
                    continue;
                }
                for (index_t k = 0; k < 2; k ++) {
                    if (i > j) 
                        root->sdoublet[k][i][j] = root->sdoublet[k][j][i];
                    else 
                        root->sdoublet[k][i][j] = 
                            root->children[0]->sdoublet[k][i][j] * root->children[0]->support_[k] +
                            root->children[1]->sdoublet[k][i][j] * root->children[1]->support_[k] + (
                                root->children[0]->ssinglet[i] * root->children[1]->ssinglet[j] +
                                root->children[1]->ssinglet[i] * root->children[0]->ssinglet[j] * (i != 0 || j != 0)) * 
                            root->children[0]->length_ * root->children[1]->length_;
                    root->mdoublet[k][i] += root->sdoublet[k][i][j];
                }
            }
            root->tdoublet[0] += root->mdoublet[0][i];
            root->tdoublet[1] += root->mdoublet[1][i];
        }
        root->tdoublet[0] = (root->tdoublet[0] - root->sdoublet[0][0][0]) / 2;
        root->tdoublet[1] = (root->tdoublet[1] - root->sdoublet[1][0][0]) / 2;
    }
}

void Tree::build_sdoublet_(Node *root) {
    if (root->parent == NULL) {
        for (index_t i = 0; i <= root->size; i ++) {
            for (index_t j = 0; j <= root->size; j ++) 
                root->sdoublet_[0][i][j] = root->sdoublet_[1][i][j] = 0;
            root->mdoublet_[0][i] = root->mdoublet_[1][i] = 0;
        }
        root->tdoublet_[0] = root->tdoublet_[1] = 0;
    }
    if (root->children.size() == 0) return ;
    for (index_t c = 0; c < 2; c ++) {
        root->children[c]->tdoublet_[0] = root->children[c]->tdoublet_[1] = 0;
        for (index_t i = 0; i <= root->size; i ++) {
            root->children[c]->mdoublet_[0][i] = root->children[c]->mdoublet_[1][i] = 0;
            for (index_t j = 0; j <= root->size; j ++) {
                if (i > 0 && i == j) {
                    root->sdoublet_[0][i][j] = 0;
                    continue;
                }
                for (index_t k = 0; k < 2; k ++) {
                    if (i > j) 
                        root->children[c]->sdoublet_[k][i][j] = root->children[c]->sdoublet_[k][j][i];
                    else 
                        root->children[c]->sdoublet_[k][i][j] = 
                            root->sdoublet_[k][i][j] * root->children[c]->support_[k] +
                            root->children[1 - c]->sdoublet[k][i][j] * root->children[1 - c]->support_[k] * root->children[c]->support_[k] + (
                                root->ssinglet_[i] * root->children[1 - c]->ssinglet[j] +
                                root->ssinglet_[j] * root->children[1 - c]->ssinglet[i] * (i != 0 || j != 0)) *
                            root->children[c]->support_[k] * root->children[1 - c]->length_;
                    root->children[c]->mdoublet_[k][i] += root->children[c]->sdoublet_[k][i][j];
                }
            }
            root->children[c]->tdoublet_[0] += root->children[c]->mdoublet_[0][i];
            root->children[c]->tdoublet_[1] += root->children[c]->mdoublet_[1][i];
        }
        root->children[c]->tdoublet_[0] = (root->children[c]->tdoublet_[0] - root->children[c]->sdoublet_[0][0][0]) / 2;
        root->children[c]->tdoublet_[1] = (root->children[c]->tdoublet_[1] - root->children[c]->sdoublet_[1][0][0]) / 2;
        build_sdoublet_(root->children[c]);
    }
}

void Tree::build_striplet(Node *root) {
    if (root->children.size() == 0) {
        for (index_t i = 0; i <= root->size; i ++) {
            for (index_t j = 0; j <= root->size; j ++) {
                root->striplet[0][i][j] = root->striplet[1][i][j] = 0;
            }
        }
    }
    else {
        for (Node *child : root->children) 
            build_striplet(child);
        for (index_t i = 1; i <= root->size; i ++) {
            for (index_t j = 0; j <= root->size; j ++) {
                if (i == j) {
                    root->striplet[0][i][j] = 0;
                    continue;
                }
                for (index_t k = 0; k < 2; k ++) {
                    root->striplet[k][i][j] = 
                        root->children[0]->striplet[k][i][j] * root->children[0]->length_ + 
                        root->children[1]->striplet[k][i][j] * root->children[1]->length_ +
                        squartet(
                            [root, k, i](index_t z) {return root->children[0]->sdoublet[k][i][z];},
                            [root, k](index_t z) {return root->children[1]->ssinglet[z];},
                            root->size, i, j
                        ) * root->children[0]->support_[k] * root->children[1]->length_ +
                        squartet(
                            [root, k](index_t z) {return root->children[0]->ssinglet[z];},
                            [root, k, i](index_t z) {return root->children[1]->sdoublet[k][i][z];},
                            root->size, i, j
                        ) * root->children[1]->support_[k] * root->children[0]->length_;
                }
            }
        }
    }
}

void Tree::build_striplet_(Node *root) {
    if (root->children.size() == 0) {
        for (index_t i = 0; i <= root->size; i ++) {
            for (index_t j = 0; j <= root->size; j ++) {
                root->striplet[0][i][j] = root->striplet[1][i][j] = 0;
            }
        }
    }
    else {
        for (Node *child : root->children) 
            build_striplet_(child);
        for (index_t i = 1; i <= root->size; i ++) {
            for (index_t j = 0; j <= root->size; j ++) {
                if (i == j) {
                    root->striplet[0][i][j] = 0;
                    continue;
                }
                for (index_t k = 0; k < 2; k ++) {
                    weight_t d0 = root->children[0]->sdoublet[k][0][0] + root->children[0]->tdoublet[k] - root->children[0]->mdoublet[k][i];
                    if (j != 0) d0 = d0 - root->children[0]->mdoublet[k][j] + root->children[0]->sdoublet[k][i][j];
                    weight_t d1 = root->children[1]->sdoublet[k][0][0] + root->children[1]->tdoublet[k] - root->children[1]->mdoublet[k][i];
                    if (j != 0) d1 = d1 - root->children[1]->mdoublet[k][j] + root->children[1]->sdoublet[k][i][j];
                    root->striplet[k][i][j] = 
                        root->children[0]->striplet[k][i][j] * root->children[0]->length_ + 
                        root->children[1]->striplet[k][i][j] * root->children[1]->length_ + 
                        root->children[0]->length_ * root->children[0]->ssinglet[i] * root->children[1]->support_[k] * d1 +
                        root->children[1]->length_ * root->children[1]->ssinglet[i] * root->children[0]->support_[k] * d0;
                }
            }
        }
    }
}

std::unordered_set<index_t> Tree::wg_edges(Node *root, Taxa &subset, weight_t ***graph) {
    if (root->children.size() == 0) {
        std::unordered_set<index_t> subtree;
        index_t index = subset.get_index(root->index);
        subtree.insert(index);
        if (subset.root_key(root->index) == 0) {
            for (index_t k = 0; k < 2; k ++) {
                for (index_t i = 0; i <= root->size; i ++) 
                    root->pdoublet[k][i] = 0;
                for (index_t i = 0; i <= root->size; i ++) 
                    root->ptriplet[k][i] = 0;
            }
            root->plength = 1;
        }
        return subtree;
    }
    else {
        Node *left = root->children[0], *right = root->children[1], *parent = root->parent;
        std::unordered_set<index_t> left_subtree = wg_edges(left, subset, graph);
        std::unordered_set<index_t> right_subtree = wg_edges(right, subset, graph);
        for (index_t i : left_subtree) {
            for (index_t j : right_subtree) {
                if (i == j) continue;
                index_t x = subset.root_key(i), y = subset.root_key(j);
                index_t i_ = subset.root_index(i), j_ = subset.root_index(j);
                weight_t s[2] = {0, 0};
                for (index_t k = 0; k < 2; k ++) {
                    if (x == 0) {
                        Node *nx = index2node[i];
                        if (y == 0) {
                            Node *ny = index2node[j];
                            s[k] += nx->ptriplet[k][0] * ny->plength * left->length_ * right->length_;
                            s[k] += squartet(
                                    [nx, k](index_t i) {return nx->pdoublet[k][i];}, 
                                    [root](index_t i) {return root->ssinglet_[i];}, 
                                    root->size, 0, 0
                                ) * ny->plength * left->support_[k] * right->length_;
                            s[k] += squartet(
                                    [nx, k](index_t i) {return nx->pdoublet[k][i];},
                                    [ny, k](index_t i) {return ny->pdoublet[k][i];},
                                    root->size, 0, 0
                                ) * left->support_[k] * right->support_[k];
                            s[k] += nx->plength * squartet(
                                    [root](index_t i) {return root->ssinglet_[i];},
                                    [ny, k](index_t i) {return ny->pdoublet[k][i];},
                                    root->size, 0, 0
                                ) * left->length_ * right->support_[k];
                            s[k] += nx->plength * ny->ptriplet[k][0] * left->length_ * right->length_;
                        }
                        else {
                            s[k] += nx->ptriplet[k][y] * right->ssinglet[y] * left->length_ * right->length_;
                            s[k] += squartet(
                                    [nx, k](index_t i) {return nx->pdoublet[k][i];},
                                    [root](index_t i) {return root->ssinglet_[i];},
                                    root->size, 0, y
                                ) * right->ssinglet[y] * left->support_[k] * right->length_;
                            s[k] += squartet(
                                    [nx, k](index_t i) {return nx->pdoublet[k][i];},
                                    [right, y, k](index_t i) {return right->sdoublet[k][y][i];},
                                    root->size, 0, y
                                ) * left->support_[k] * right->support_[k];
                            s[k] += nx->plength * squartet(
                                    [root](index_t i) {return root->ssinglet_[i];},
                                    [right, y, k](index_t i) {return right->sdoublet[k][y][i];},
                                    root->size, 0, y
                                ) * left->length_ * right->support_[k];
                            s[k] += nx->plength * right->striplet[k][y][0] * left->length_ * right->length_;
                        }
                    }
                    else {
                        if (y == 0) {
                            Node *ny = index2node[j];
                            s[k] += left->striplet[k][x][0] * ny->plength * left->length_ * right->length_;
                            s[k] += squartet(
                                    [left, x, k](index_t i) {return left->sdoublet[k][x][i];},
                                    [root](index_t i) {return root->ssinglet_[i];},
                                    root->size, x, 0
                                ) * ny->plength * left->support_[k] * right->length_;
                            s[k] += squartet(
                                    [left, x, k](index_t i) {return left->sdoublet[k][x][i];},
                                    [ny, k](index_t i) {return ny->pdoublet[k][i];},
                                    root->size, x, 0
                                ) * left->support_[k] * right->support_[k];
                            s[k] += left->ssinglet[x] * squartet(
                                    [root](index_t i) {return root->ssinglet_[i];},
                                    [ny, k](index_t i) {return ny->pdoublet[k][i];},
                                    root->size, x, 0
                                ) * left->length_ * right->support_[k];
                            s[k] += left->ssinglet[x] * ny->ptriplet[k][x] * left->length_ * right->length_;
                        }
                        else {
                            s[k] += left->striplet[k][x][y] * right->ssinglet[y] * left->length_ * right->length_;
                            s[k] += squartet(
                                    [left, x, k](index_t i) {return left->sdoublet[k][x][i];},
                                    [root](index_t i) {return root->ssinglet_[i];},
                                    root->size, x, y
                                ) * right->ssinglet[y] * left->support_[k] * right->length_;
                            s[k] += squartet(
                                    [left, x, k](index_t i) {return left->sdoublet[k][x][i];},
                                    [right, y, k](index_t i) {return right->sdoublet[k][y][i];},
                                    root->size, x, y
                                ) * left->support_[k] * right->support_[k];
                            s[k] += left->ssinglet[x] * squartet(
                                    [root](index_t i) {return root->ssinglet_[i];},
                                    [right, y, k](index_t i) {return right->sdoublet[k][y][i];},
                                    root->size, x, y
                                ) * left->length_ * right->support_[k];
                            s[k] += left->ssinglet[x] * right->striplet[k][y][x] * left->length_ * right->length_;
                        }
                    }
                }
                graph[0][i_][j_] += s[1] - s[0];
                graph[0][j_][i_] += s[1] - s[0]; 
            }
        }
        std::unordered_set<index_t> subtree;
        for (index_t i : left_subtree) {
            subtree.insert(i);
            if (subset.root_key(i) == 0) {
                Node *nx = index2node[i];
                for (index_t c = 0; c < 2; c ++) {
                    for (index_t k = 0; k <= root->size; k ++) {
                        nx->ptriplet[c][k] = nx->ptriplet[c][k] * left->length_ +
                            squartet(
                                [nx, c](index_t i) {return nx->pdoublet[c][i];},
                                [right](index_t i) {return right->ssinglet[i];},
                                root->size, 0, k
                            ) * left->support_[c] * right->length_;
                    }
                    for (index_t k = 0; k <= root->size; k ++) {
                        nx->pdoublet[c][k] = nx->pdoublet[c][k] * left->support_[c] +
                            nx->plength * right->ssinglet[k] * left->length_ * right->length_;
                    }
                }
                nx->plength *= left->length_;
            }
        }
        for (index_t j : right_subtree) {
            subtree.insert(j);
            if (subset.root_key(j) == 0) {
                Node *ny = index2node[j];
                for (index_t c = 0; c < 2; c ++) {
                    for (index_t k = 0; k <= root->size; k ++) {
                        ny->ptriplet[c][k] = ny->ptriplet[c][k] * right->length_ +
                            squartet(
                                [left](index_t i) {return left->ssinglet[i];},
                                [ny, c](index_t i) {return ny->pdoublet[c][i];},
                                root->size, 0, k
                            ) * left->length_ * right->support_[c];
                    }
                    for (index_t k = 0; k <= root->size; k ++) {
                        ny->pdoublet[c][k] = ny->pdoublet[c][k] * right->support_[c] +
                            left->ssinglet[k] * ny->plength * left->length_ * right->length_;
                    }
                }
                ny->plength *= right->length_;
            }
        }
        if (DEBUG_MODE) test_pxlet(root, subtree, subset);
        return subtree;
    }
}

std::unordered_set<index_t> Tree::wb_edges(Node *root, Taxa &subset, weight_t ***graph) {
    if (root->children.size() == 0) {
        std::unordered_set<index_t> subtree;
        index_t index = subset.get_index(root->index);
        subtree.insert(index);
        if (subset.root_key(root->index) == 0) {
            for (index_t k = 0; k < 2; k ++) {
                for (index_t i = 0; i <= root->size; i ++) 
                    root->ptriplet[k][i] = 0;
            }
            root->plength = 1;
        }
        return subtree;
    }
    else {
        Node *left = root->children[0], *right = root->children[1], *parent = root->parent;
        std::unordered_set<index_t> left_subtree = wb_edges(left, subset, graph);
        std::unordered_set<index_t> right_subtree = wb_edges(right, subset, graph);
        for (index_t i : left_subtree) {
            for (index_t j : right_subtree) {
                if (i == j) continue;
                index_t x = subset.root_key(i), y = subset.root_key(j);
                index_t i_ = subset.root_index(i), j_ = subset.root_index(j);
                weight_t s[2] = {0, 0};
                for (index_t k = 0; k < 2; k ++) {
                    if (x == 0) {
                        Node *nx = index2node[i];
                        if (y == 0) {
                            Node *ny = index2node[j];
                            s[k] += nx->ptriplet[k][0] * left->length_ * ny->plength * right->length_;
                            s[k] += nx->plength * left->length_ * ny->plength * right->length_ *
                                (root->sdoublet_[k][0][0] + root->tdoublet_[k]);
                            s[k] += nx->plength * left->length_ * ny->ptriplet[k][0] * right->length_;
                        }
                        else {
                            s[k] += nx->ptriplet[k][y] * left->length_ * right->ssinglet[y] * right->length_;
                            s[k] += nx->plength * left->length_ * right->ssinglet[y] * right->length_ *
                                (root->sdoublet_[k][0][0] + root->tdoublet_[k] - root->mdoublet_[k][y]);
                            s[k] += nx->plength * left->length_ * right->striplet[k][y][0] * right->length_;
                        }
                    }
                    else {
                        if (y == 0) {
                            Node *ny = index2node[j];
                            s[k] += left->striplet[k][x][0] * left->length_ * ny->plength * right->length_;
                            s[k] += left->ssinglet[x] * left->length_ * ny->plength * right->length_ *
                                (root->sdoublet_[k][0][0] + root->tdoublet_[k] - root->mdoublet_[k][x]);
                            s[k] += left->ssinglet[x] * left->length_ * ny->ptriplet[k][x] * right->length_;
                        }
                        else {
                            s[k] += left->striplet[k][x][y] * left->length_ * right->ssinglet[y] * right->length_;
                            s[k] += left->ssinglet[x] * left->length_ * right->ssinglet[y] * right->length_ * (
                                root->sdoublet_[k][0][0] + root->tdoublet_[k] -
                                root->mdoublet_[k][x] - root->mdoublet_[k][y] + root->sdoublet_[k][x][y]);
                            s[k] += left->ssinglet[x] * left->length_ * right->striplet[k][y][x] * right->length_;
                        }
                    }
                }
                graph[1][i_][j_] += s[1] - s[0];
                graph[1][j_][i_] += s[1] - s[0]; 
            }
        }
        std::unordered_set<index_t> subtree;
        for (index_t i : left_subtree) {
            subtree.insert(i);
            if (subset.root_key(i) == 0) {
                Node *nx = index2node[i];
                for (index_t c = 0; c < 2; c ++) {
                    for (index_t k = 0; k <= root->size; k ++) {
                        weight_t d = right->sdoublet[c][0][0] + right->tdoublet[c];
                        if (k != 0) d -= right->mdoublet[c][k];
                        nx->ptriplet[c][k] = 
                            nx->ptriplet[c][k] * left->length_ + 
                            nx->plength * left->length_ * d * right->support_[c];
                    }
                }
                nx->plength *= left->length_;
            }
        }
        for (index_t j : right_subtree) {
            subtree.insert(j);
            if (subset.root_key(j) == 0) {
                Node *ny = index2node[j];
                for (index_t c = 0; c < 2; c ++) {
                    for (index_t k = 0; k <= root->size; k ++) {
                        weight_t d = left->sdoublet[c][0][0] + left->tdoublet[c];
                        if (k != 0) d-= left->mdoublet[c][k];
                        ny->ptriplet[c][k] = 
                            ny->ptriplet[c][k] * right->length_ + 
                            ny->plength * right->length_ * d * left->support_[c];
                    }
                }
                ny->plength *= right->length_;
            }
        }
        if (DEBUG_MODE) test_pxlet_(root, subtree, subset);
        return subtree;
    }
}

void Tree::test_pxlet(Node *root, std::unordered_set<index_t> &subtree, Taxa &subset) {
    for (index_t k = 0; k < 2; k ++) {
        for (index_t i : subtree) {
            if (subset.root_key(i) == 0) {
                Node *nx = index2node[i];
                Node *a = nx;
                weight_t t = 1;
                while (a != root) {t *= a->length_; a = a->parent;}
                weight_t err = fabs(t - nx->plength);
                if (err > ERROR_BOUND) {
                    std::cout << "p1 " << (double)t << ' ' << (double)err << std::endl;
                }
                assert(err <= ERROR_BOUND);
                std::vector<Node *> leaves;
                get_leaves(root, &leaves);
                for (index_t j = 0; j <= root->size; j ++) {
                    weight_t s = 0;
                    for (Node* ny : leaves) {
                        if (subset.root_key(ny->index) != j) continue;
                        if (ny == nx) continue;
                        Node *a = nx, *b = ny;
                        weight_t t = 1;
                        while (a->depth > b->depth) {t *= a->length_; a = a->parent;}
                        while (b->depth > a->depth) {t *= b->length_; b = b->parent;}
                        while (a != b) {t *= a->length_ * b->length_; a = a->parent; b = b->parent;}
                        while (a != root) {t *= a->support_[k]; a = a->parent;}
                        s += t * subset.root_weight(ny->index);
                    }
                    weight_t err = fabs(s - nx->pdoublet[k][j]);
                    if (err > ERROR_BOUND) {
                        std::cout << "p2 " << (double)t << ' ' << (double)err << std::endl;
                    }
                    assert(err <= ERROR_BOUND);
                }
                for (index_t j = 0; j <= root->size; j ++) {
                    weight_t s = 0;
                    Node *x = nx;
                    for (Node* y : leaves) {
                        for (Node *z : leaves) {
                            index_t yk = subset.root_key(y->index);
                            index_t zk = subset.root_key(z->index);
                            if (y == x || z == x) continue;
                            if (zk != 0 && (yk == zk)) continue;
                            if (zk == 0 && (y == z)) continue;
                            if (j != 0 && (yk == j || zk == j)) continue;
                            Node *a = y, *b = z;
                            while (a->depth > b->depth) {a = a->parent;}
                            while (b->depth > a->depth) {b = b->parent;}
                            while (a != b) {a = a->parent; b = b->parent;}
                            Node *yz = a;
                            a = x; b = y;
                            while (a->depth > b->depth) {a = a->parent;}
                            while (b->depth > a->depth) {b = b->parent;}
                            while (a != b) {a = a->parent; b = b->parent;}
                            Node *xy = a;
                            a = x; b = z;
                            while (a->depth > b->depth) {a = a->parent;}
                            while (b->depth > a->depth) {b = b->parent;}
                            while (a != b) {a = a->parent; b = b->parent;}
                            Node *xz = a;
                            if (yz->depth > xy->depth && yz->depth > xz->depth) continue;
                            weight_t t = 1;
                            if (xy->depth > xz->depth) {
                                a = x;
                                while (a != xy) {t *= a->length_; a = a->parent;}
                                a = y;
                                while (a != xy) {t *= a->length_; a = a->parent;}
                                a = z;
                                while (a != xz) {t *= a->length_; a = a->parent;}
                                a = xy;
                                while (a != xz) {t *= a->support_[k]; a = a->parent;}
                                a = xz;
                                while (a != root) {t *= a->length_; a = a->parent;}
                            }
                            else {
                                a = x;
                                while (a != xz) {t *= a->length_; a = a->parent;}
                                a = z;
                                while (a != xz) {t *= a->length_; a = a->parent;}
                                a = y;
                                while (a != xy) {t *= a->length_; a = a->parent;}
                                a = xz;
                                while (a != xy) {t *= a->support_[k]; a = a->parent;}
                                a = xy;
                                while (a != root) {t *= a->length_; a = a->parent;}
                            }
                            s += t * subset.root_weight(y->index) * subset.root_weight(z->index);
                        }
                    }
                    s /= 2;
                    weight_t err = fabs(s - nx->ptriplet[k][j]);
                    if (err > ERROR_BOUND) {
                        std::cout << "p3 " << (double)s << ' ' << (double)err << std::endl;
                    }
                    assert(err <= ERROR_BOUND);
                }
            }
        }
    }
}

void Tree::test_pxlet_(Node *root, std::unordered_set<index_t> &subtree, Taxa &subset) {
    for (index_t k = 0; k < 2; k ++) {
        for (index_t i : subtree) {
            if (subset.root_key(i) == 0) {
                Node *nx = index2node[i];
                Node *a = nx;
                weight_t t = 1;
                while (a != root) {t *= a->length_; a = a->parent;}
                weight_t err = fabs(t - nx->plength);
                if (err > ERROR_BOUND) {
                    std::cout << "p1_ " << (double)t << ' ' << (double)err << std::endl;
                }
                assert(err <= ERROR_BOUND);
                std::vector<Node *> leaves;
                get_leaves(root, &leaves);
                for (index_t j = 0; j <= root->size; j ++) {
                    weight_t s = 0;
                    Node *x = nx;
                    for (Node* y : leaves) {
                        for (Node *z : leaves) {
                            index_t yk = subset.root_key(y->index);
                            index_t zk = subset.root_key(z->index);
                            if (y == x || z == x) continue;
                            if (zk != 0 && (yk == zk)) continue;
                            if (zk == 0 && (y == z)) continue;
                            if (j != 0 && (yk == j || zk == j)) continue;
                            Node *a = y, *b = z;
                            while (a->depth > b->depth) {a = a->parent;}
                            while (b->depth > a->depth) {b = b->parent;}
                            while (a != b) {a = a->parent; b = b->parent;}
                            Node *yz = a;
                            a = x; b = y;
                            while (a->depth > b->depth) {a = a->parent;}
                            while (b->depth > a->depth) {b = b->parent;}
                            while (a != b) {a = a->parent; b = b->parent;}
                            Node *xy = a;
                            a = x; b = z;
                            while (a->depth > b->depth) {a = a->parent;}
                            while (b->depth > a->depth) {b = b->parent;}
                            while (a != b) {a = a->parent; b = b->parent;}
                            Node *xz = a;
                            if (yz->depth > xy->depth && yz->depth > xz->depth) {
                                weight_t t = 1;
                                Node *a = y, *b = z;
                                while (a->depth > b->depth) {t *= a->length_; a = a->parent;}
                                while (b->depth > a->depth) {t *= b->length_; b = b->parent;}
                                while (a != b) {t *= a->length_; t *= b->length_; a = a->parent; b = b->parent;}
                                b = x;
                                while (a->depth > b->depth) {t *= a->support_[k]; a = a->parent;}
                                while (b->depth > a->depth) {t *= b->length_; b = b->parent;}
                                while (a != b) {t *= a->support_[k]; t *= b->length_; a = a->parent; b = b->parent;}
                                while (a != root) {t *= a->length_; a = a->parent;}
                                s += t * subset.root_weight(y->index)  * subset.root_weight(z->index);
                            }
                        }
                    }
                    s /= 2;
                    weight_t err = fabs(s - nx->ptriplet[k][j]);
                    if (err > ERROR_BOUND) {
                        std::cout << "p3_ " << (double)s << ' ' << (double)err << ' ' << (double)nx->ptriplet[k][j] << std::endl;
                    }
                    assert(err <= ERROR_BOUND);
                }
            }
        }
    }
}

weight_t ***Tree::build_wgraph(Taxa &subset) {
    if (DEBUG_MODE) get_depth(root, 0);
    build_wstates(root, subset);
    build_ssinglet(root, subset);
    if (DEBUG_MODE) test_ssinglet(root, subset);
    build_ssinglet_(root);
    if (DEBUG_MODE) test_ssinglet_(root, subset);
    build_sdoublet(root);
    if (DEBUG_MODE) test_sdoublet(root, subset);
    build_striplet(root);
    if (DEBUG_MODE) test_striplet(root, subset);
    weight_t ***graph = new weight_t**[2];
    graph[0] = Matrix::new_mat(subset.size());
    graph[1] = Matrix::new_mat(subset.size());
    wg_edges(root, subset, graph);
    build_sdoublet_(root);
    if (DEBUG_MODE) test_sdoublet_(root, subset);
    build_striplet_(root);
    if (DEBUG_MODE) test_striplet_(root, subset);
    wb_edges(root, subset, graph);
    if (DEBUG_MODE) test_graph(root, subset, graph);
    clear_wstates(root);
    return graph;
}

void Tree::test_graph(Node *root, Taxa &subset, weight_t ***graph) {
    std::unordered_map<quartet_t, weight_t> quartets;
    get_wquartets(&quartets);
    std::unordered_map<index_t, index_t> index2index;
    for (index_t i = 0; i < subset.size(); i ++) {
        index2index[subset.root_at(i)] = i;
    }
    weight_t ***bf_graph = new weight_t**[2];
    bf_graph[0] = Matrix::new_mat(subset.size());
    bf_graph[1] = Matrix::new_mat(subset.size());
    for (auto elem : quartets) {
        index_t *indices = split(elem.first);
        weight_t I = 1;
        for (index_t i = 0; i < 4; i ++) {
            I *= subset.root_weight(indices[i]);
            indices[i] = subset.get_index(indices[i]);
        }
        for (index_t i = 0; i < 4; i ++) 
            for (index_t j = i + 1; j < 4; j ++) 
                if (indices[i] == indices[j]) I = 0;
        index_t a = index2index[indices[0]], b = index2index[indices[1]], c = index2index[indices[2]], d = index2index[indices[3]];
        weight_t w = elem.second * I;
        bf_graph[1][a][b] += w; bf_graph[1][c][d] += w; bf_graph[1][b][a] += w; bf_graph[1][d][c] += w;
        bf_graph[0][a][c] += w; bf_graph[0][a][d] += w; bf_graph[0][b][c] += w; bf_graph[0][b][d] += w;
        bf_graph[0][c][a] += w; bf_graph[0][d][a] += w; bf_graph[0][c][b] += w; bf_graph[0][d][b] += w;
        delete [] indices;
    }
    /*
    std::cout << subset.size() << std::endl;
    std::cout << Matrix::display_mat(bf_graph[0], subset.size()) << std::endl;
    std::cout << Matrix::display_mat(graph[0], subset.size()) << std::endl;
    std::cout << std::endl;
    */
    for (index_t i = 0; i < subset.size(); i ++) {
        for (index_t j = 0; j < subset.size(); j ++) {
            if (i == j) continue;
            weight_t err = fabs(graph[0][i][j] - bf_graph[0][i][j]);
            if (err > ERROR_BOUND) {
                std::cout << "g0 " << (double)bf_graph[0][i][j] << ' ' << (double)err << ' ' << (double)graph[0][i][j] << std::endl;
            }
            assert(err <= ERROR_BOUND);
            err = fabs(graph[1][i][j] - bf_graph[1][i][j]);
            if (err > ERROR_BOUND) {
                std::cout << "g1 " << (double)bf_graph[1][i][j] << ' ' << (double)err << ' ' << (double)graph[1][i][j] << std::endl;
                std::cout << subset.singleton_taxa() << std::endl;
                std::cout << i << ' ' << j << std::endl;
            }
            assert(err <= ERROR_BOUND);
        }
    }
    Matrix::delete_mat(bf_graph[0], subset.size());
    Matrix::delete_mat(bf_graph[1], subset.size());
    delete [] bf_graph;
}

void Tree::test_ssinglet(Node *root, Taxa &subset) {
    for (Node *child : root->children) 
        test_ssinglet(child, subset);
    std::vector<Node *> leaves;
    get_leaves(root, &leaves);
    for (index_t i = 0; i <= root->size; i ++) {
        weight_t s = 0;
        for (Node *leaf : leaves) {
            Node *p = leaf;
            weight_t t = 0;
            while (p != root) {
                t += p->length;
                p = p->parent;
            }
            index_t index = subset.root_key(leaf->index);
            if (index == i)
                s += subset.root_weight(leaf->index) * exp(- t);
        }
        weight_t eps = fabs(root->ssinglet[i] - s);
        if (eps > ERROR_BOUND) 
            std::cout << "s1 " << (double)eps << ' ' << (double)root->ssinglet[i] << ' ' << (double)s << std::endl;
        assert(eps <= ERROR_BOUND);
    }
}

void Tree::test_ssinglet_(Node *root, Taxa &subset) {
    for (Node *child : root->children) 
        test_ssinglet_(child, subset);
    std::vector<Node *> leaves;
    get_leaves(root, &leaves);
    std::unordered_set<Node *> leaf_set;
    for (Node *leaf : leaves) leaf_set.insert(leaf);
    leaves.clear();
    get_leaves(this->root, &leaves);
    for (index_t i = 0; i <= root->size; i ++) {
        weight_t s = 0;
        for (Node *leaf : leaves) {
            if (leaf_set.find(leaf) != leaf_set.end()) continue;
            Node *a = leaf, *b = root;
            weight_t t = 0;
            while (a->depth > b->depth) {t += a->length; a = a->parent;}
            while (b->depth > a->depth) {t += b->length; b = b->parent;}
            while (a != b) {t += a->length + b->length; a = a->parent; b = b->parent;}
            index_t index = subset.root_key(leaf->index);
            if (index == i)
                s += subset.root_weight(leaf->index) * exp(- t);
        }
        weight_t eps = fabs(root->ssinglet_[i] - s);
        if (eps > ERROR_BOUND) 
            std::cout << "s1_ " << (double)eps << std::endl;
        assert(eps <= ERROR_BOUND);
    }
}

void Tree::test_sdoublet(Node *root, Taxa &subset) {
    for (Node *child : root->children) 
        test_sdoublet(child, subset);
    for (index_t k = 0; k < 2; k ++) {
        std::vector<Node *> leaves;
        get_leaves(root, &leaves);
        for (index_t i = 0; i <= root->size; i ++) {
            for (index_t j = 0; j <= root->size; j ++) {
                if (i > 0 && i == j) continue;
                weight_t s = 0;
                for (Node *x : leaves) {
                    for (Node *y : leaves) {
                        if (subset.root_key(x->index) != i) continue;
                        if (subset.root_key(y->index) != j) continue;
                        if (x == y) continue;
                        Node *a = x, *b = y;
                        weight_t t = 1;
                        while (a->depth > b->depth) {t *= a->length_; a = a->parent;}
                        while (b->depth > a->depth) {t *= b->length_; b = b->parent;}
                        while (a != b) {t *= a->length_ * b->length_; a = a->parent; b = b->parent;}
                        while (a != root) {t *= a->support_[k]; a = a->parent;}
                        s += t * subset.root_weight(x->index) * subset.root_weight(y->index);
                    }
                }
                if (i == 0 && j == 0) s /= 2;
                weight_t err = fabs(s - root->sdoublet[k][i][j]);
                if (err > ERROR_BOUND) {
                    std::cout << "s2 " << (double)s << ' ' << (double)err << ' ' << (double)root->sdoublet[k][i][j]<< ' ' << i << ' ' << j << std::endl;
                }
                assert(err <= ERROR_BOUND);
            }
        }
    }
}

void Tree::test_sdoublet_(Node *root, Taxa &subset) {
    std::vector<Node *> leaves;
    get_leaves(root, &leaves);
    std::unordered_set<Node *> leaf_set;
    for (Node *leaf : leaves) leaf_set.insert(leaf);
    leaves.clear();
    get_leaves(this->root, &leaves);
    for (index_t k = 0; k < 2; k ++) {
        for (index_t i = 0; i <= root->size; i ++) {
            for (index_t j = 0; j <= root->size; j ++) {
                if (i > 0 && i == j) continue;
                weight_t s = 0;
                for (Node *x : leaves) {
                    if (leaf_set.find(x) != leaf_set.end()) continue;
                    for (Node *y : leaves) {
                        if (leaf_set.find(y) != leaf_set.end()) continue;
                        if (subset.root_key(x->index) != i) continue;
                        if (subset.root_key(y->index) != j) continue;
                        if (x == y) continue;
                        Node *a = x, *b = y;
                        weight_t t = 1;
                        while (a->depth > b->depth) {t *= a->length_; a = a->parent;}
                        while (b->depth > a->depth) {t *= b->length_; b = b->parent;}
                        while (a != b) {t *= a->length_ * b->length_; a = a->parent; b = b->parent;}
                        Node *xy = a;
                        a = x; b = root;
                        while (a->depth > b->depth) {a = a->parent;}
                        while (b->depth > a->depth) {b = b->parent;}
                        while (a != b) {a = a->parent; b = b->parent;}
                        Node *xr = a;
                        a = y; b = root;
                        while (a->depth > b->depth) {a = a->parent;}
                        while (b->depth > a->depth) {b = b->parent;}
                        while (a != b) {a = a->parent; b = b->parent;}
                        Node *yr = a;
                        if (xr == yr) {
                            Node *a = xy, *b = root;
                            while (a->depth > b->depth) {t *= a->support_[k]; a = a->parent;}
                            while (b->depth > a->depth) {t *= b->support_[k]; b = b->parent;}
                            while (a != b) {t *= a->support_[k]; t *= b->support_[k]; a = a->parent; b = b->parent;}
                        }
                        if (xy == xr) {
                            Node *a = root;
                            while (a != yr) {t *= a->support_[k]; a = a->parent;}
                        }
                        if (xy == yr) {
                            Node *a = root;
                            while (a != xr) {t *= a->support_[k]; a = a->parent;}
                        }
                        s += t * subset.root_weight(x->index) * subset.root_weight(y->index);
                    }
                }
                if (i == 0 && j == 0) s /= 2;
                weight_t err = fabs(s - root->sdoublet_[k][i][j]);
                if (err > ERROR_BOUND) {
                    std::cout << "s2_ " << (double)s << ' ' << (double)err << ' ' << (double)root->sdoublet_[k][i][j] << ' ' << i << ' ' << j << std::endl;
                }
                assert(err <= ERROR_BOUND);
            }
        }
    }
    for (Node *child : root->children) 
        test_sdoublet_(child, subset);
}

void Tree::test_striplet(Node *root, Taxa &subset) {
    for (Node *child : root->children) 
        test_striplet(child, subset);
    for (index_t k = 0; k < 2; k ++) {
        std::vector<Node *> leaves;
        get_leaves(root, &leaves);
        for (index_t i = 1; i <= root->size; i ++) {
            for (index_t j = 0; j <= root->size; j ++) {
                if (i == j) continue;
                weight_t s = 0;
                for (Node *x : leaves) {
                    if (subset.root_key(x->index) != i) continue;
                    for (Node *y : leaves) {
                        for (Node *z : leaves) {
                            index_t yk = subset.root_key(y->index);
                            index_t zk = subset.root_key(z->index);
                            if (zk != 0 && (yk == zk)) continue;
                            if (zk == 0 && (y == z)) continue;
                            if (yk == i || zk == i) continue;
                            if (j != 0 && (yk == j || zk == j)) continue;
                            Node *a = y, *b = z;
                            while (a->depth > b->depth) {a = a->parent;}
                            while (b->depth > a->depth) {b = b->parent;}
                            while (a != b) {a = a->parent; b = b->parent;}
                            Node *yz = a;
                            a = x; b = y;
                            while (a->depth > b->depth) {a = a->parent;}
                            while (b->depth > a->depth) {b = b->parent;}
                            while (a != b) {a = a->parent; b = b->parent;}
                            Node *xy = a;
                            a = x; b = z;
                            while (a->depth > b->depth) {a = a->parent;}
                            while (b->depth > a->depth) {b = b->parent;}
                            while (a != b) {a = a->parent; b = b->parent;}
                            Node *xz = a;
                            if (yz->depth > xy->depth && yz->depth > xz->depth) continue;
                            weight_t t = 1;
                            if (xy->depth > xz->depth) {
                                a = x;
                                while (a != xy) {t *= a->length_; a = a->parent;}
                                a = y;
                                while (a != xy) {t *= a->length_; a = a->parent;}
                                a = z;
                                while (a != xz) {t *= a->length_; a = a->parent;}
                                a = xy;
                                while (a != xz) {t *= a->support_[k]; a = a->parent;}
                                a = xz;
                                while (a != root) {t *= a->length_; a = a->parent;}
                            }
                            else {
                                a = x;
                                while (a != xz) {t *= a->length_; a = a->parent;}
                                a = z;
                                while (a != xz) {t *= a->length_; a = a->parent;}
                                a = y;
                                while (a != xy) {t *= a->length_; a = a->parent;}
                                a = xz;
                                while (a != xy) {t *= a->support_[k]; a = a->parent;}
                                a = xy;
                                while (a != root) {t *= a->length_; a = a->parent;}
                            }
                            s += t * subset.root_weight(x->index) * subset.root_weight(y->index)  * subset.root_weight(z->index);
                        }
                    }
                }
                s /= 2;
                weight_t err = fabs(s - root->striplet[k][i][j]);
                if (err > ERROR_BOUND) {
                    std::cout << "s3 " << (double)s << ' ' << (double)err << std::endl;
                }
                assert(err <= ERROR_BOUND);
            }
        }
    }
}

void Tree::test_striplet_(Node *root, Taxa &subset) {
    for (Node *child : root->children) 
        test_striplet_(child, subset);
    for (index_t k = 0; k < 2; k ++) {
        std::vector<Node *> leaves;
        get_leaves(root, &leaves);
        for (index_t i = 1; i <= root->size; i ++) {
            for (index_t j = 0; j <= root->size; j ++) {
                if (i == j) continue;
                weight_t s = 0;
                for (Node *x : leaves) {
                    if (subset.root_key(x->index) != i) continue;
                    for (Node *y : leaves) {
                        for (Node *z : leaves) {
                            index_t yk = subset.root_key(y->index);
                            index_t zk = subset.root_key(z->index);
                            if (zk != 0 && (yk == zk)) continue;
                            if (zk == 0 && (y == z)) continue;
                            if (yk == i || zk == i) continue;
                            if (j != 0 && (yk == j || zk == j)) continue;
                            Node *a = y, *b = z;
                            while (a->depth > b->depth) {a = a->parent;}
                            while (b->depth > a->depth) {b = b->parent;}
                            while (a != b) {a = a->parent; b = b->parent;}
                            Node *yz = a;
                            a = x; b = y;
                            while (a->depth > b->depth) {a = a->parent;}
                            while (b->depth > a->depth) {b = b->parent;}
                            while (a != b) {a = a->parent; b = b->parent;}
                            Node *xy = a;
                            a = x; b = z;
                            while (a->depth > b->depth) {a = a->parent;}
                            while (b->depth > a->depth) {b = b->parent;}
                            while (a != b) {a = a->parent; b = b->parent;}
                            Node *xz = a;
                            if (yz->depth > xy->depth && yz->depth > xz->depth) {
                                weight_t t = 1;
                                Node *a = y, *b = z;
                                while (a->depth > b->depth) {t *= a->length_; a = a->parent;}
                                while (b->depth > a->depth) {t *= b->length_; b = b->parent;}
                                while (a != b) {t *= a->length_; t *= b->length_; a = a->parent; b = b->parent;}
                                b = x;
                                while (a->depth > b->depth) {t *= a->support_[k]; a = a->parent;}
                                while (b->depth > a->depth) {t *= b->length_; b = b->parent;}
                                while (a != b) {t *= a->support_[k]; t *= b->length_; a = a->parent; b = b->parent;}
                                while (a != root) {t *= a->length_; a = a->parent;}
                                s += t * subset.root_weight(x->index) * subset.root_weight(y->index)  * subset.root_weight(z->index);
                            }
                        }
                    }
                }
                s /= 2;
                weight_t err = fabs(s - root->striplet[k][i][j]);
                if (err > ERROR_BOUND) {
                    std::cout << "s3 " << (double)s << ' ' << (double)err << ' ' << i << ' ' << j << std::endl;
                }
                assert(err <= ERROR_BOUND);
            }
        }
    }
}
