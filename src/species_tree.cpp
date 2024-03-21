#include "tree.hpp"
#include "graph.hpp"
#include "toms743.hpp"

SpeciesTree::SpeciesTree(std::string stree_file, Dict *dict) {
    this->dict = dict;
    this->artifinyms = dict->max_size();
    this->root = NULL;
    //this->leaf_for_rooting = NULL;
    //this->qfreq_mode = "";

    std::unordered_map<std::string, std::string> dummy;
    std::string newick;
    std::size_t tracker = dict->size();

    // Read species tree
    std::cout << "Reading species tree" << std::endl;

    std::ifstream fin(stree_file);
    if (fin.fail()) {
        std::cout << "\nERROR: Unable to open species tree file " << stree_file << std::endl;
        exit(1);
    }
    std::getline(fin, newick);
    // TODO: Check validity of newick string
    if (newick.find(";") == std::string::npos)
        std::cout << "\nERROR: Not a valid newick string" << std::endl;
    root = build_tree(newick, dummy);
    fin.close();

    if (dict->size() > tracker) {
        std::cout << "\nERROR: Species tree has more taxa than input!" << std::endl;
        exit(1);
    }

    if (this->size() < 3) {
        std::cout << "\nERROR: Species tree has less than 4 taxa!" << std::endl;
        exit(1);
    }

    std::cout << "Found" << std::endl;
    std::cout << "    " << size() << " taxa" << std::endl;

    // TODO: Add suppress uniforcations!

    // Refined polytomies
    int total = refine_tree(root);
    if (total > 0)
        std::cout << "    " << total << " polytomies" << std::endl;

}

SpeciesTree::SpeciesTree(std::vector<Tree *> &input, Dict *dict, std::string mode, unsigned long int iter_limit, std::string output_file) {
    this->dict = dict;
    this->artifinyms = dict->max_size();
    this->mode = mode;
    this->iter_limit = iter_limit;
    this->root = NULL;
    //this->leaf_for_rooting = NULL;
    //this->qfreq_mode = mode[3];

    Taxa subset(dict, mode);
    switch (mode[1]) {
        case '0': {
            // Fast execution mode
            std::cout << "Constructing species tree" << std::endl;
            std::cout << "At most " << subset.size() * 2 - 5 << " subproblems:" << std::endl;
            root = construct_stree(input, subset, -1, 0);
            std::cout << std::setw(5) << count[0] << " subproblems computed." << std::endl;
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
            std::ofstream fout(output_file + "_quartets.txt");
            if (fout.fail()) {
                std::cout << "\nERROR: Unable to write to file " << output_file << "_quartets.txt" << std::endl;
                exit(1);
            }

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

            for (auto elem : quartets) {
                index_t *indices = split(elem.first);
                index_t a = indices[0], b = indices[1], c = indices[2], d = indices[3];
                fout << "(("
                     << dict->index2label(a) << ',' << dict->index2label(b) 
                     << "),("
                     << dict->index2label(c) << ',' << dict->index2label(d)
                     << ")); "
                     << std::fixed << std::setprecision(16) << (double)elem.second << std::endl;
                delete [] indices;
            }

            fout.close();

            break;
        }
        case '3': {
            // Compute good and bad edges, then exit
            std::ofstream good_edges_txt(output_file + "_good_edges.txt");
            if (good_edges_txt.fail()) {
                std::cout << "\nERROR: Unable to write to file " << output_file << "_good_edges.txt" << std::endl;
                exit(1);
            }

            std::ofstream bad_edges_txt(output_file + "_bad_edges.txt");
            if (bad_edges_txt.fail()) {
                std::cout << "\nERROR: Unable to write to file " << output_file << "_bad_edges.txt" << std::endl;
                exit(1);
            }

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

void SpeciesTree::get_qfreq_around_branch(Node *root, std::vector<Tree *> &input, std::string &qfreq_mode) {
    // Check tree is binary
    if (root->children.size() != 2 && root->children.size() != 0) {
        std::cout << "ERROR: Cannot compute quartet support in non-binary tree" << std::endl;
        exit(1);
    }

    // Initalize quartet frequencies to 0
    root->f[0] = 0;
    root->f[1] = 0;
    root->f[2] = 0;

    // Check if there are quartets around the branch
    bool ok = true;
    if (root->parent == NULL) {
        // At root
        ok = false;
    }
    else if (root->is_leaf()) {
        // At leaf
        // TODO: handle multi-individuals
        ok = false;
    } else if (root->parent->parent == NULL) {
        // At node incident to root
        if (root->get_sibling()->is_leaf()) {
            // Rooted at leaf so skip
            ok = false;
        }
        else if (root != this->root->children[0]) {
            // Not rooted at leaf so pick one branch
            ok = false;
        }
    }

    // Compute quartet frequencies        
    if (ok) {
        std::vector<weight_t> localf0, localf1, localf2, localqs;
        std::unordered_map<index_t, index_t> quad;
        std::vector<Node *> x, y, z, w;

        Node *mysib = root->get_sibling();
        get_leaves(this->root, &w);           // get all leaves in tree
        get_leaves(root->children[0], &x);    // get leaves below left child
        get_leaves(root->children[1], &y);    // get leaves below right child
        if (root->parent->parent == NULL) {
            // At root's left child
            get_leaves(mysib->children[0], &z);
        } else {
            get_leaves(mysib, &z);  // get leaves below sibling
        }

        // Let w' = w / x,y,z
        // First quartet topology x,y|z,w'
        for (Node *leaf : w) quad[leaf->index] = 4;
        for (Node *leaf : x) quad[leaf->index] = 1; // sibling
        for (Node *leaf : y) quad[leaf->index] = 2; // sibling
        for (Node *leaf : z) quad[leaf->index] = 3;
        for (Tree *t : input)
            localf0.push_back(t->get_qfreq(quad) / 2);

        // Second quartet topology x,z|y,w'
        for (Node *leaf : w) quad[leaf->index] = 4;
        for (Node *leaf : x) quad[leaf->index] = 1; // sibling
        for (Node *leaf : y) quad[leaf->index] = 3;
        for (Node *leaf : z) quad[leaf->index] = 2; // sibling
        for (Tree *t : input)
            localf1.push_back(t->get_qfreq(quad) / 2);

        // Third quartet topology x,w'|y,z
        for (Node *leaf : w) quad[leaf->index] = 2; // sibling
        for (Node *leaf : x) quad[leaf->index] = 1; // sibling
        for (Node *leaf : y) quad[leaf->index] = 3;
        for (Node *leaf : z) quad[leaf->index] = 4;
        for (Tree *t : input)
            localf2.push_back(t->get_qfreq(quad) / 2);

        if ((qfreq_mode == "n") || (qfreq_mode == "f")) {
            // Processing for normal (unweighted) quartet frequencies

            // Compute the total quartet weight per gene tree because
            // it can be impacted by missing data and polytomies
            for (Tree *t : input)
                localqs.push_back(t->get_qcount(quad));

            // Combine quartet frequencies across gene trees
            for (index_t i = 0; i < input.size(); i++) {
                // Normalize so each gene tree votes once (but can split its vote)
                if (localqs[i] != 0.0) {
                    localf0[i] /= localqs[i];  // x,y|z,w'
                    localf1[i] /= localqs[i];  // x,z|y,w'
                    localf2[i] /= localqs[i];  // x,w'|y,z
                    // frequency of 1,2,3,4 is (localqs[i] - (localf0[i] + localf1[i] + localf2[i])) / localqs[i];
                }

                // Sum votes
                root->f[0] += localf0[i];
                root->f[1] += localf1[i];
                root->f[2] += localf2[i];
            }
        }
        else {
            // Processing for weighted quartet frequencies (similar to weighted ASTRAL)
            weight_t localtotal, sumsq = 0.0, gtcount = 0.0;

            for (std::size_t i = 0; i < input.size(); i++) {
                localtotal = localf0[i] + localf1[i] + localf2[i];
                sumsq += localtotal * localtotal;
                if ((localtotal) > 1e-9) gtcount++;
                root->f[0] += localf0[i];
                root->f[1] += localf1[i];
                root->f[2] += localf2[i];
            }

            if (sumsq > 0) {
                root->f[0] *= sqrt(gtcount / sumsq);
                root->f[1] *= sqrt(gtcount / sumsq);
                root->f[2] *= sqrt(gtcount / sumsq);
            }
        }
    }

    for (Node *child : root->children) {
        get_qfreq_around_branch(child, input, qfreq_mode);
    }
}

void SpeciesTree::write_pcs_table(std::vector<Tree *> &input, std::vector<std::size_t> &positions, std::string &qfreq_mode, std::ostream &os) {
    if (pcs_node == NULL) {
        std::cout << "ERROR: Did not find node in species tree for computing PCS" << std::endl;
        exit(1);
    }

    // Check if there are quartets around the branch
    bool ok = true;
    if (pcs_node->parent == NULL) {
        ok = false;
    }
    else if (pcs_node->is_leaf()) {
        // TODO: handle multi-individuals
        ok = false;
    } else if (pcs_node->parent->parent == NULL) {
        if (pcs_node->get_sibling()->is_leaf()) {
            ok = false;
        }
    }

    if (!ok) {
        std::cout << "ERROR: Cannot compute PCS bceause no quartets around flagged branch" << std::endl;
        exit(1);
    }

    // Compute quartet frequencies    
    std::unordered_map<index_t, index_t> quad;
    std::vector<Node *> x, y, z, w;
    std::unordered_set<Node *> wprime;
    Node *mysib = pcs_node->get_sibling();

    get_leaves(pcs_node->children[0], &x);    // get leaves below left child
    get_leaves(pcs_node->children[1], &y);    // get leaves below right child
    if (pcs_node->parent->parent == NULL) {
        // At root's child
        get_leaves(mysib->children[0], &z);
        get_leaves(mysib->children[1], &w);
        get_leaf_set(mysib->children[1], &wprime);
    } else {
        get_leaves(mysib, &z);             // get leaves below sibling
        get_leaves(this->root, &w);        // get all leaves in tree and take complement
        get_leaf_set(this->root, &wprime);   
    }
    for (Node *leaf : x) wprime.erase(leaf);
    for (Node *leaf : y) wprime.erase(leaf);
    for (Node *leaf : z) wprime.erase(leaf);

    std::vector<weight_t> localf0, localf1, localf2, localqs;
    weight_t localtotal;

    // Let w' = w / x,y,z
    // First quartet topology x,y|z,w'
    for (Node *leaf : w) quad[leaf->index] = 4;
    for (Node *leaf : x) quad[leaf->index] = 1;  // sibling
    for (Node *leaf : y) quad[leaf->index] = 2;  // sibling
    for (Node *leaf : z) quad[leaf->index] = 3;
    for (Tree *t : input)
        localf0.push_back(t->get_qfreq(quad) / 2);

    // Second quartet topology x,z|y,w'
    for (Node *leaf : w) quad[leaf->index] = 4;
    for (Node *leaf : x) quad[leaf->index] = 1;  // sibling
    for (Node *leaf : y) quad[leaf->index] = 3;
    for (Node *leaf : z) quad[leaf->index] = 2;  // sibling
    for (Tree *t : input)
        localf1.push_back(t->get_qfreq(quad) / 2);

    // Third quartet topology x,w'|y,z
    for (Node *leaf : w) quad[leaf->index] = 2;  // sibling
    for (Node *leaf : x) quad[leaf->index] = 1;  // sibling
    for (Node *leaf : y) quad[leaf->index] = 3;
    for (Node *leaf : z) quad[leaf->index] = 4;
    for (Tree *t : input)
        localf2.push_back(t->get_qfreq(quad) / 2);

    os << "# x = ";
    for (Node* leaf : x) os << dict->index2label(leaf->index) << ",";
    os << "\b\n# y = ";
    for (Node* leaf : y) os << dict->index2label(leaf->index) << ",";
    os << "\b\n# z = ";
    for (Node* leaf : z) os << dict->index2label(leaf->index) << ",";
    os << "\b\n# w = ";
    for (Node* leaf : wprime) os << dict->index2label(leaf->index) << ",";
    os << "\b\n";

    if ((qfreq_mode == "n") || (qfreq_mode == "f")) {
        // Processing for normal (unweighted) quartet frequencies

        os << "id\tposition\tf_xy|zw\tf_xz|yw\tf_xw|yz\tf_xyzw\ttotalf\n";

        // Check the fourth topology by figuring out
        // quartet count per gene tree and subtracting 
        for (Tree *t : input)
            localqs.push_back(t->get_qcount(quad));

        // Write to CSV
        for (index_t i = 0; i < input.size(); i++) {
            localtotal = localf0[i] + localf1[i] + localf2[i];


            //std::cout << "here " << localtotal << std::endl;
           
            if (localtotal > 0) {
                os << i << "\t"
                   << positions[i] << "\t"
                   << localf0[i] << "\t"
                   << localf1[i] << "\t"
                   << localf2[i] << "\t"
                   << localqs[i] - localtotal << "\t"
                   << localqs[i] << "\n";
            }
        }
    }
    else {
        // Processing for weighted quartet frequencies (similar to weighted ASTRAL)
        weight_t f0, f1, f2, sumsq = 0.0, gtcount = 0.0;

        os << "id\tposition\tf_xy|zw\tf_xz|yw\tf_xw|yz\ttotalf\n";

        // Initialize
        pcs_node->f[0] = 0;
        pcs_node->f[1] = 0;
        pcs_node->f[2] = 0;

        // Complete calculation for normalization
        for (std::size_t i = 0; i < input.size(); i++) {
            localtotal = localf0[i] + localf1[i] + localf2[i];
            sumsq += localtotal * localtotal;
            if ((localtotal) > 1e-9) gtcount++;
            pcs_node->f[0] += localf0[i];
            pcs_node->f[1] += localf1[i];
            pcs_node->f[2] += localf2[i];
        }

        // Write to CSV
        for (index_t i = 0; i < input.size(); i++) {
            f0 = localf0[i] * sqrt(gtcount / sumsq);
            f1 = localf1[i] * sqrt(gtcount / sumsq);
            f2 = localf2[i] * sqrt(gtcount / sumsq);
            localtotal = f0 + f1 + f2;

            if (localtotal > 0) {
                os << i << "\t"
                   << positions[i] << "\t"
                   << f0 << "\t"
                   << f1 << "\t"
                   << f2 << "\t"
                   << localtotal << "\n";
            }
        }
    }
}

std::string compute_branch_length_msc_gt(weight_t q1) {
    if (q1 < (1.0 / 3.0))
        return "0";

    if (q1 == 1.0)
        return "9";  // same as MP-EST
    
    weight_t x = -1.0 * log(1.5 * (1.0 - q1));

    return std::to_string(x);
}

std::string compute_branch_length_msc_bp(weight_t q1) {
    if (q1 < (1.0 / 3.0))
        return "0";

    if (q1 == 1.0)
        return "9";  // same as MP-EST

    weight_t y, x;
    int outcome;

    y = ((2.0 / 3.0) / (1.0 - q1)) - 1.0;
    x = toms743::wapr(y, 0, outcome, 0);

    if (outcome != 1)
        return std::to_string(x);
    else
        return "nan";
}


void SpeciesTree::annotate(std::vector<Tree *> input, std::string &qfreq_mode) {
    get_qfreq_around_branch(root, input, qfreq_mode);
}

std::string SpeciesTree::to_string_annotated(std::string brln_mode) {
    return display_tree_annotated(root, brln_mode);
}

std::string SpeciesTree::display_tree_annotated(Node *root, std::string brln_mode) {
    if (root->is_leaf()) {
        // TODO: update for multi-individual setting
        return dict->index2label(root->index);
    }

    std::string s = "(";
    for (Node * node : root->children)
        s += display_tree_annotated(node, brln_mode) + ",";
    s[s.size() - 1] = ')';

    // Check if there are quartets around the branch (incoming edge to root)
    if (root->parent == NULL) {
        return s + ";";
    }
    else if (root->parent->parent == NULL) {
        if (root->get_sibling()->is_leaf()) {
            // Rooted at leaf so don't annotate
            return s;
        }
        else if (root != this->root->children[0]) {
            // Not rooted at leaf so pick one branch to annotate
            return s;
        }
    }

    // If quartets around branch (incoming edge to root), write info
    weight_t f1, f2, f3, en, q1 = 0, q2 = 0, q3 = 0;
    f1 = root->f[0];
    f2 = root->f[1];
    f3 = root->f[2];
    en = f1 + f2 + f3;
    if (en > 0) {
        q1 = f1 / en;
        q2 = f2 / en;
        q3 = f3 / en;
    }

    // Compute branch length from q1
    std::string brlen = "";
    if (brln_mode == "g")
        brlen = ":" + compute_branch_length_msc_gt(q1);
    else if (brln_mode == "b")
        brlen = ":" + compute_branch_length_msc_bp(q1);

    // Return string for node with branch info
    return s + "\'[f1=" + std::to_string(f1) +
                 ";f2=" + std::to_string(f2) +
                 ";f3=" + std::to_string(f3) +
                 ";q1=" + std::to_string(q1) +
                 ";q2=" + std::to_string(q2) +
                 ";q3=" + std::to_string(q3) +
                 ";EN=" + std::to_string(en) +
                 "]\'" + brlen;
}

void SpeciesTree::write_support_table_row(Node *root, std::ostream &os, std::string brln_mode) {
    if (root->children.size() != 2 && root->children.size() != 0) {
        std::cout << "ERROR: Cannot write table for non-binary tree" << std::endl;
        exit(1);
    }

    // Check if there are quartets around the branch
    bool ok = true;
    if (root->parent == NULL) {
        // At root
        ok = false;
    }
    else if (root->is_leaf()) {
        // At leaf
        // TODO: handle multi-individuals
        ok = false;
    } else if (root->parent->parent == NULL) {
        // At node incident to root
        if (root->get_sibling()->is_leaf()) {
            // Rooted at leaf so skip
            ok = false;
        } else if (root != this->root->children[0]) {
            // Not rooted at leaf so pick one branch
            ok = false;
        }
    }

    if (ok) {
        // Get taxa around branch
        std::unordered_set<Node *> x, y, z, w;
        Node *mysib = root->get_sibling();

        get_leaf_set(root->children[0], &x);    // get leaves below left child
        get_leaf_set(root->children[1], &y);    // get leaves below right child
        if (root->parent->parent == NULL) {
            // At root's left child
            get_leaf_set(mysib->children[0], &z);
            get_leaf_set(mysib->children[1], &w);
        } else {
            get_leaf_set(mysib, &z);            // get leaves below sibling
            get_leaf_set(this->root, &w);       // get all leaves in tree and take complement
            for (Node *leaf : x) w.erase(leaf);
            for (Node *leaf : y) w.erase(leaf);
            for (Node *leaf : z) w.erase(leaf);
        }

        // Get quartet support
        weight_t f1, f2, f3, en, q1 = 0, q2 = 0, q3 = 0;
        f1 = root->f[0];
        f2 = root->f[1];
        f3 = root->f[2];
        en = f1 + f2 + f3;
        if (en > 0) {
            q1 = f1 / en;
            q2 = f2 / en;
            q3 = f3 / en;
        }

        // Get branch length
        std::string brlen = "NA";
        if (brln_mode == "g")
            brlen = compute_branch_length_msc_gt(q1);
        else if (brln_mode == "b")
            brlen = compute_branch_length_msc_bp(q1);

        // Write row for t1 = x,y|z,w'
        os << "N" << root->index << "\tt1\t{";
        for (Node* leaf : z) os << dict->index2label(leaf->index) << ",";
        os << "\b}|{";
        for (Node* leaf : x) os << dict->index2label(leaf->index) << ",";
        os << "\b}#{";
        for (Node* leaf : y) os << dict->index2label(leaf->index) << ",";
        os << "\b}|{";
        for (Node* leaf : w) os << dict->index2label(leaf->index) << ",";
        os << "\b}\tNA\t" << f1 << "\t" << en << "\t" << q1 << "\t" << brlen << std::endl;

        // Write row for t2 = x,z|y,w'
        os << "N" << root->index << "\tt2\t{";
        for (Node* leaf : y) os << dict->index2label(leaf->index) << ",";
        os << "\b}|{";
        for (Node* leaf : x) os << dict->index2label(leaf->index) << ",";
        os << "\b}#{";
        for (Node* leaf : z) os << dict->index2label(leaf->index) << ",";
        os << "\b}|{";
        for (Node* leaf : w) os << dict->index2label(leaf->index) << ",";
        os << "\b}\tNA\t" << f2 << "\t" << en << "\t" << q2 << "\tNA" << std::endl;

        // Write row for t3 = x,w'|y,z
        os << "N" << root->index << "\tt3\t{";
        for (Node* leaf : y) os << dict->index2label(leaf->index) << ",";
        os << "\b}|{";
        for (Node* leaf : x) os << dict->index2label(leaf->index) << ",";
        os << "\b}#{";
        for (Node* leaf : w) os << dict->index2label(leaf->index) << ",";
        os << "\b}|{";
        for (Node* leaf : z) os << dict->index2label(leaf->index) << ",";
        os << "\b}\tNA\t" << f3 << "\t" << en << "\t" << q3 << "\tNA" << std::endl;
    }

    for (Node *child : root->children) {
        write_support_table_row(child, os, brln_mode);
    }
}

void SpeciesTree::write_support_table(std::ostream &os, std::string brln_mode) {
    os << "NODE\tTREE\tQUARTET\tSKIP\tQFREQ\nEN\tQSUPP\tBRLEN\n";
    write_support_table_row(root, os, brln_mode);
}

void SpeciesTree::root_at_clade(std::unordered_set<std::string> &clade_label_set) {
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
            std::cout << "    Instead rooting at " << backup_label << std::endl;
            node = find_node(backup_index);
        }
    }

    if (node == NULL) {
        std::cout << "    Unable to root tree" << std::endl;
        return;
    }
    reroot_on_edge_above_node(node);
}

/*void SpeciesTree::put_back_root() {
    if (leaf_for_rooting != NULL)
        reroot_on_edge_above_node(leaf_for_rooting);
}*/



