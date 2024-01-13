#include "tree.hpp"
#include "graph.hpp"

SpeciesTree::SpeciesTree(std::string stree_file, Dict *dict) {
    this->dict = dict;
    this->artifinyms = dict->max_size();
    leaf_for_rooting = NULL;

    std::ifstream fin(stree_file);
    if (fin.fail()) {
        std::cout << std::endl << "ERROR: Unable to read file " << stree_file << std::endl;
        exit(1);
    }

    std::string newick;
    index_t i = 0;
    while (std::getline(fin, newick)) {
        if (newick.find(";") == std::string::npos) break;
        root = build_tree(newick);
    }
    fin.close();

    if (this->size() < 3) {
        std::cout << "ERROR: Score tree has less than 4 leaves!" << std::endl;
        exit(1);
    }
    std::cout << "Score tree has " << size() << " leaves." << std::endl;


    // Resolve polytomies
    int total = resolve_tree(root);
    if (total > 0)
        std::cout << "Resolved " << total << " polytomies in " << stree_file << std::endl;

    // Re-position root
    Node *sibling, *node_for_rooting;
    for (Node* child: root->children) {
        if (child->children.size() == 0)
            leaf_for_rooting = child;
        else 
            sibling = child;
    }

    if (leaf_for_rooting != NULL) {
        // If rooted at leaf, balance around root so get_freq function works
        for (Node* child: sibling->children) {
            if (child->children.size() >= 2)
                node_for_rooting = child;
        }
        reroot_on_edge_above_node(node_for_rooting);
    }
}

SpeciesTree::SpeciesTree(std::vector<Tree *> &input, Dict *dict, std::string mode, unsigned long int iter_limit) {
    this->dict = dict;
    this->artifinyms = dict->max_size();
    this->mode = mode;
    this->iter_limit = iter_limit;
    leaf_for_rooting = NULL;

    Taxa subset(dict, mode);
    switch (mode[1]) {
        case '0': {
            // Fast execution mode
            std::cout << "At most " << subset.size() * 2 - 5 << " subproblems.\n";
            root = construct_stree(input, subset, -1, 0);
            std::cout << std::setw(5) << count[0] << " subproblems computed.\n";
            break;
        }
        case '1': {
            // Brute force execution mode
            std::unordered_map<quartet_t, weight_t> quartets;
            if (mode[3] == 'f') {
                // Use unweighted code
                for (Tree * tree: input) tree->get_quartets(&quartets);
            } else if (mode[3] == 'n' || mode[3] == 's') {
                // Use weighted support only code
                for (Tree * tree: input) tree->get_wquartets_(&quartets);
            } else {
                // Use weighted hybrid code; also used for length only
                for (Tree * tree: input) tree->get_wquartets(&quartets);
            }
            root = construct_stree(quartets, subset, -1, 0);
            break;
        }
        case '2': {
            // Compute weighted quartets, then exit
            std::unordered_map<quartet_t, weight_t> quartets;
            if (mode[3] == 'f') {
                // Use unweighted code
                for (Tree * tree: input) tree->get_quartets(&quartets);
            } else if (mode[3] == 'n' || mode[3] == 's') {
                // Use weighted support only code
                for (Tree * tree: input) tree->get_wquartets_(&quartets);
            } else {
                // Use weighted hybrid code; also used for length only
                for (Tree * tree: input) tree->get_wquartets(&quartets);
            }
            // std::cout << to_string(quartets);
            // std::cout << quartets.size() << std::endl;
            // if (quartets_txt.is_open()) {
                for (auto elem : quartets) {
                    index_t *indices = split(elem.first);
                    index_t a = indices[0], b = indices[1], c = indices[2], d = indices[3];
                    quartets_txt
                        << "(("
                        << dict->index2label(a) << ',' << dict->index2label(b) 
                        << "),("
                        << dict->index2label(c) << ',' << dict->index2label(d)
                        << ")); "
                        << std::fixed << std::setprecision(16) << (double)elem.second << std::endl;
                    delete [] indices;
                }
                quartets_txt.close();
            //}
            break;
        }
        case '3': {
            // Compute good and bad edges, then exit
            Graph *g = new Graph(input, subset, mode.substr(3, 1));

            g->write_good_edges(dict);
            good_edges_txt.close();

            g->write_bad_edges(dict);
            bad_edges_txt.close();
            break;
        }
        default: {
            break;
        }
    }
    // std::cout << artifinyms << std::endl;
}

SpeciesTree::~SpeciesTree() {
    
}

index_t SpeciesTree::artifinym() {
    return -- artifinyms;
}

Node *SpeciesTree::construct_stree(std::vector<Tree *> &input, Taxa &subset, index_t parent_pid, index_t depth) {
    if (count[0] % 10 == 0) std::cout << std::setw(5) << count[0] << " subproblems computed.\n";
    //std::cout << "Computing subproblem:" << std::setw(5) << count[0] << std::flush;
    //std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
    index_t size = subset.size(), pid = count[0] ++;
    if (verbose > "0") subproblem_csv << pid << ',' << parent_pid << ',' << depth;
    if (verbose > "0") subproblem_csv << ',' << subset.size() << ',' << subset.artificial_taxa();
    if (verbose > "0") subproblem_csv << ',' << "\"" + subset.to_list() + "\"";
    Node *root;
    if (size < 4) {
        if (size == 1) {
            root = new Node(subset.root_at(0));
        }
        else if (size == 2) {
            root = new Node(pseudonym());
            root->children.push_back(new Node(subset.root_at(0)));
            root->children.push_back(new Node(subset.root_at(1)));
            root->children[0]->parent = root->children[1]->parent = root;
        }
        else {
            root = new Node(pseudonym());
            root->children.push_back(new Node(pseudonym()));
            Node *left = root->children[0];
            left->children.push_back(new Node(subset.root_at(0)));
            left->children.push_back(new Node(subset.root_at(1)));
            left->children[0]->parent = left->children[1]->parent = left;
            root->children.push_back(new Node(subset.root_at(2)));
            root->children[0]->parent = root->children[1]->parent = root;
        }
        if (verbose > "0") subproblem_csv << std::endl;
    }
    else {
        /*
        for (auto tree : input) {
            std::unordered_map<index_t, index_t> valid = tree->get_indices();
            subset.weight_update(valid);
            weight_t ***subgraph = tree->build_wgraph(subset);
            Matrix::delete_mat(subgraph[0], size);
            Matrix::delete_mat(subgraph[1], size);
            delete [] subgraph;
        }
        */
        Graph *g = new Graph(input, subset, mode.substr(3, 1));
        if (verbose > "0") subproblem_csv << std::endl;
        std::vector<index_t> A, B;
        weight_t max = g->get_cut(&A, &B, iter_limit);
        if (max < 0) {
            root = new Node(pseudonym());
            for (index_t i = 0; i < subset.size(); i ++) {
                Node *child = new Node(subset.root_at(i));
                root->children.push_back(child);
                child->parent = root;
            }
        }
        else {
            Taxa subsetA(subset), subsetB(subset);
            index_t artificial = artifinym();
            subsetA.struct_update(A, artificial);
            subsetB.struct_update(B, artificial);
            root = new Node(pseudonym());
            root->children.push_back(reroot_stree(construct_stree(input, subsetA, pid, depth + 1), artificial));
            root->children.push_back(reroot_stree(construct_stree(input, subsetB, pid, depth + 1), artificial));
            root->children[0]->parent = root->children[1]->parent = root;
        }
        delete g;
    }
    // std::cout << display_tree(root) << std::endl;
    return root;
}

Node *SpeciesTree::construct_stree(std::unordered_map<quartet_t, weight_t> &quartets, Taxa &subset, index_t parent_pid, index_t depth) {
    index_t size = subset.size(), pid = count[0] ++;
    if (verbose > "0") subproblem_csv << pid << ',' << parent_pid << ',' << depth;
    if (verbose > "0") subproblem_csv << ',' << subset.size() << ',' << subset.artificial_taxa();
    if (verbose > "0") subproblem_csv << ',' << "\"" + subset.to_list() + "\"";
    Node *root;
    if (size < 4) {
        if (size == 1) {
            root = new Node(subset.root_at(0));
        }
        else if (size == 2) {
            root = new Node(pseudonym());
            root->children.push_back(new Node(subset.root_at(0)));
            root->children.push_back(new Node(subset.root_at(1)));
            root->children[0]->parent = root->children[1]->parent = root;
        }
        else {
            root = new Node(pseudonym());
            root->children.push_back(new Node(pseudonym()));
            Node *left = root->children[0];
            left->children.push_back(new Node(subset.root_at(0)));
            left->children.push_back(new Node(subset.root_at(1)));
            left->children[0]->parent = left->children[1]->parent = left;
            root->children.push_back(new Node(subset.root_at(2)));
            root->children[0]->parent = root->children[1]->parent = root;
        }
        if (verbose > "0") subproblem_csv << std::endl;
    }
    else {
        if (verbose > "0") subproblem_csv << std::endl;
        // std::cout << to_string(quartets) << std::endl;
        Graph *g = new Graph(quartets, subset);
        std::vector<index_t> A, B;
        weight_t max = g->get_cut(&A, &B, iter_limit);
        if (max < 0) {
            root = new Node(pseudonym());
            for (index_t i = 0; i < subset.size(); i ++) {
                Node *child = new Node(subset.root_at(i));
                root->children.push_back(child);
                child->parent = root;
            }
        }
        else {
            std::unordered_set<index_t> setA(A.begin(), A.end()), setB(B.begin(), B.end());
            Taxa subsetA(subset), subsetB(subset);
            index_t artificial = artifinym();
            subsetA.struct_update(B, artificial);
            subsetB.struct_update(A, artificial);
            root = new Node(pseudonym());
            std::unordered_map<quartet_t, weight_t> quartetsA, quartetsB;
            for (auto elem : quartets) {
                index_t *indices = split(elem.first);
                index_t count = 0;
                for (index_t i = 0; i < 4; i ++) {
                    if (setA.find(indices[i]) != setA.end()) 
                        count ++;
                }
                switch (count) {
                    case 0: {
                        if (quartetsB.find(elem.first) == quartetsB.end()) 
                            quartetsB[elem.first] = 0;
                        quartetsB[elem.first] += elem.second;
                        break;
                    }
                    case 1: {
                        for (index_t i = 0; i < 4; i ++) {
                            if (setA.find(indices[i]) != setA.end()) {
                                indices[i] = artificial;
                            }
                        }
                        quartet_t temp = join(indices);
                        if (quartetsB.find(temp) == quartetsB.end()) 
                            quartetsB[temp] = 0;
                        quartetsB[temp] += elem.second / A.size();
                        break;
                    }
                    case 2: {
                        break;
                    }
                    case 3: {
                        for (index_t i = 0; i < 4; i ++) {
                            if (setB.find(indices[i]) != setB.end()) {
                                indices[i] = artificial;
                            }
                        }
                        quartet_t temp = join(indices);
                        if (quartetsA.find(temp) == quartetsA.end()) 
                            quartetsA[temp] = 0;
                        quartetsA[temp] += elem.second / B.size();
                        break;
                    }
                    case 4: {
                        if (quartetsA.find(elem.first) == quartetsA.end()) 
                            quartetsA[elem.first] = 0;
                        quartetsA[elem.first] += elem.second;
                        break;
                    }
                }
                delete [] indices;
            }
            root->children.push_back(reroot_stree(construct_stree(quartetsA, subsetA, pid, depth + 1), artificial));
            root->children.push_back(reroot_stree(construct_stree(quartetsB, subsetB, pid, depth + 1), artificial));
            root->children[0]->parent = root->children[1]->parent = root;
            delete g;
        }
    }
    // std::cout << display_tree(root) << std::endl;
    return root;
}

Node *SpeciesTree::reroot(Node *root, std::unordered_set<index_t> &visited) {
    std::vector<Node *> child;
    if (root->parent != NULL && visited.find(root->parent->index) == visited.end()) {
        visited.insert(root->parent->index);
        child.push_back(reroot(root->parent, visited));
    }
    for (Node *ch : root->children) {
        if (ch != NULL && visited.find(ch->index) == visited.end()) {
            visited.insert(ch->index);
            child.push_back(reroot(ch, visited));
        }
    }
    Node *new_root;
    if (child.size() >= 2) {
        new_root = new Node(pseudonym());
        visited.insert(new_root->index);
        for (Node *ch : child) {
            new_root->children.push_back(ch);
            ch->parent = new_root;
        }
        return new_root;
    }
    else if (child.size() == 1) {
        new_root = child[0];
    }
    else {
        new_root = new Node(root->index);
    }
    // std::cout << '>' << display_tree(new_root) << std::endl;
    return new_root;
}

Node *SpeciesTree::reroot_stree(Node *root, index_t artificial) {
    Node *new_root = artificial2node(root, artificial);
    std::unordered_set<index_t> visited;
    visited.insert(new_root->index);
    Node *new_tree = reroot(new_root, visited);
    delete root;
    return new_tree;
}

Node *SpeciesTree::artificial2node(Node *root, index_t artificial) {
    if (root->children.size() == 0) {
        if (root->index == artificial) 
            return root;
        return NULL;
    }
    else {
        for (Node *child : root->children) {
            Node *temp = artificial2node(child, artificial);
            if (temp != NULL) return temp;
        }
        return NULL;
    }
}

void SpeciesTree::get_freq(Node *root, std::vector<Tree *> input) {
    if (root->children.size() != 0 && root->parent != NULL) {
        std::vector<Node *> x, y, z, w;
        get_leaves(root->children[0], &x);
        get_leaves(root->children[1], &y);

        if (root->parent->parent != NULL) {
            if (root->parent->children[0] == root) 
                get_leaves(root->parent->children[1], &z);
            else 
                get_leaves(root->parent->children[0], &z);
        }
        else {
            if (root->parent->children[0] == root) 
                get_leaves(root->parent->children[1]->children[0], &z);
            else 
                get_leaves(root->parent->children[0]->children[0], &z);
        }
        get_leaves(this->root, &w);

        std::vector<weight_t> localf0, localf1, localf2, localf3;
        weight_t globaltotal, localtotal, tmp;
        std::unordered_map<index_t, index_t> quad;

        root->wq[0] = 0;
        root->wq[1] = 0;
        root->wq[2] = 0;

        // First topology
        for (Node *leaf : w) quad[leaf->index] = 4;
        for (Node *leaf : x) quad[leaf->index] = 1;
        for (Node *leaf : y) quad[leaf->index] = 2;
        for (Node *leaf : z) quad[leaf->index] = 3;
        for (Tree *t : input) {
            tmp = t->get_freq(quad) / 2;;
            localf0.push_back(tmp);
            root->wq[0] += tmp;
        }

        // Second topology
        for (Node *leaf : w) quad[leaf->index] = 4;
        for (Node *leaf : x) quad[leaf->index] = 1;
        for (Node *leaf : y) quad[leaf->index] = 3;
        for (Node *leaf : z) quad[leaf->index] = 2;
        for (Tree *t : input) {
            tmp = t->get_freq(quad) / 2;
            localf1.push_back(tmp);
            root->wq[1] += tmp;
        }

        // Third topology
        for (Node *leaf : w) quad[leaf->index] = 2;
        for (Node *leaf : x) quad[leaf->index] = 1;
        for (Node *leaf : y) quad[leaf->index] = 3;
        for (Node *leaf : z) quad[leaf->index] = 4;
        for (Tree *t : input) {
            tmp = t->get_freq(quad) / 2;
            localf2.push_back(tmp);
            root->wq[2] += tmp;
        }

        // Check the fourth topology by figuring out
        // per gene tree QC and subtracting 
        for (Tree *t : input) {
            tmp = t->get_qc(quad);
            localf3.push_back(tmp);
        }

        // handle weighted astral values
        globaltotal = root->wq[0] + root->wq[1] + root->wq[2];
        if (globaltotal != 0.0) {
            root->wq[0] /= globaltotal;
            root->wq[1] /= globaltotal;
            root->wq[2] /= globaltotal;
        }

        // handle regular astral values
        for (index_t i = 0; i < input.size(); i++) {
            localtotal = localf0[i] + localf1[i] + localf2[i];
            if (localtotal != 0.0) {
                localf0[i] /= localf3[i];
                localf1[i] /= localf3[i];
                localf2[i] /= localf3[i];
            }
        }

        // Store frequencies
        root->f[0] = 0;
        root->f[1] = 0;
        root->f[2] = 0;
        for (index_t i = 0; i < input.size(); i++) {
            root->f[0] += localf0[i];
            root->f[1] += localf1[i];
            root->f[2] += localf2[i];
        }


    }
    for (Node *child : root->children) {
        get_freq(child, input);
    }
}

void SpeciesTree::annotate(std::vector<Tree *> input) {
    display_tree_index(root);
    get_freq(root, input);
    //return display_tree_annotated(root);
}

std::string SpeciesTree::to_string_annotated() {
    return display_tree_annotated(root);
}

std::string SpeciesTree::display_tree_annotated(Node *root) {
    if (root->children.size() == 0) 
        return dict->index2label(root->index);

    std::string s = "(";
    for (Node * node : root->children)
        s += display_tree_annotated(node) + ",";
    s[s.size() - 1] = ')';

    if (root->parent == NULL)
        return s + ";";


    weight_t f1, f2, f3, q1, q2, q3, en;
    std::string brlen;

    f1 = root->f[0];
    f2 = root->f[1];
    f3 = root->f[2];
    en = f1 + f2 + f3;
    q1 = f1 / en;
    q2 = f2 / en;
    q3 = f3 / en;

    if (q1 < q2 || q1 < q3)
        brlen = "-1";
    else if (q1 == 1.0)
        brlen = "40";
    else
        brlen = std::to_string(-1.0 * log(1.5 * (1 - q1)));

    return s + "\'[f1=" + std::to_string(f1) +
             ";f2=" + std::to_string(f2) +
             ";f3=" + std::to_string(f3) +
             ";q1=" + std::to_string(q1) +
             ";q2=" + std::to_string(q2) +
             ";q3=" + std::to_string(q3) +
             ";EN=" + std::to_string(en) +
             ";wq1=" + std::to_string(root->wq[0]) +
             ";wq2=" + std::to_string(root->wq[1]) +
             ";wq3=" + std::to_string(root->wq[2]) + 
             "]\':"  + brlen;

    //return s + "\'[f1=" + std::to_string(root->f[0]) + 
    //              ";f2=" + std::to_string(root->f[1]) + 
    //              ";f3=" + std::to_string(root->f[2]) + "]\'";
}

void SpeciesTree::root_at_clade(std::unordered_set<std::string> &clade_label_set) {
    std::cout << "Rooting tree" << std::endl;

    std::unordered_set<index_t> clade_index_set;
    std::string backup_label;
    Node *node;
    index_t index, backup_index;
    index_t tracker = dict->size();

    for (auto label : clade_label_set) {
        index = dict->label2index(label);
        if (dict->size() > tracker) {
            std::cout << "    " << label << " does not exist in input trees!" << std::endl;
            tracker++;
        }
        else {
            clade_index_set.insert(index);
            backup_index = index;
            backup_label = label;
        }
    }

    if (clade_index_set.size() == 0) {
        std::cout << "    None of the root taxa do not exist in tree!" << std::endl;
        return;
    }
    else if (clade_index_set.size() == 1) {
        node = find_node(backup_index);
    } else {
        node = find_node_for_split(clade_index_set);
        if (node == NULL) {
            std::cout << "    Root taxa do not form a clade in the output tree!" << std::endl;
            std::cout << "    Indead rooting at " << backup_label << std::endl;
            node = find_node(backup_index);
        }
    }

    reroot_on_edge_above_node(node);
}

void SpeciesTree::put_back_root() {
    if (leaf_for_rooting != NULL)
        reroot_on_edge_above_node(leaf_for_rooting);
}
