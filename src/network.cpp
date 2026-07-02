#include "network.hpp"
#include<queue>

size_t Network::hybrid_pseudonyms = 0;

double stod__(std::string s) {
    //std::cout << 's' << s << std::endl;
    return std::stod(s);
}

std::string to_string_p(double value, int p) {
    std::ostringstream out;
    out.precision(p);
    out << std::fixed << value;
    return std::move(out).str();
}

Network::Network(std::string &newick) {
    root = build_network(newick, label2hybrid);
    for (auto elem : label2hybrid) {
        elem.second->hybrid = label2node[elem.first];
        label2node[elem.first]->hybrid = elem.second;
    }
}

Network::Network(Network &network) {
    std::unordered_map<NetworkNode *, std::pair<std::size_t, std::size_t>> visited;
    depth_first_search(network.root, NULL, 0, visited);
    root = build_tree_of_blob(network.root);
    root = simplify_tree(root);
}

Network::Network(Tree *tree, Dict *dict) {

    std::unordered_set<Node *> candidates_node;
    std::unordered_set<index_t> candidates;
    tree->get_leaf_set(tree->get_root(), &candidates_node);
    
    for (const Node* n : candidates_node) {
        candidates.insert(n->index);
    }


    index_t candidate_index = -1;

    Node *hybrid_tree_root = tree->get_root();
    
    std::cout << "Checking whether the current rooting is valid ..." << std::endl;
    //check whether the root is null 
    if (hybrid_tree_root == nullptr || hybrid_tree_root == NULL) {
        std::cout << "The input tree is empty; will not write the output" << std::endl;
        root = NULL;
        return;
    }
    
    if (check_root(tree->get_root(), candidates)) {
        std::cout << "Cureent rooting is valid " << std::endl;
    } else {

        if (candidates.size() > 0) {
            candidate_index = *(candidates.begin());
            std::cout << "Current rooting is not valid but it is valid semi-directed network; " << std::endl;
            std::cout << "Perform rerooting on the edge above taxon " << dict->index2label(candidate_index) << std::endl;
        } else {
                std::cout << "Current network is not a valid semi-directed network; will not write the output" << std::endl;
                root=NULL; 
        }

        hybrid_tree_root = reroot_hybrid_tree_on_edge_above_node(tree->index2node[candidate_index], tree);
    
    }

    

    root = build_network(hybrid_tree_root, dict, label2hybrid);

}


Network::~Network() {
    delete root;
}

std::string Network::to_string() {
    return display_network(root) + ";";
}

std::string Network::to_string_basic() {
    return display_network_basic(root) + ";";
}

bool Network::check_root(Node *root, std::unordered_set<index_t> &candidates) {
    if (root->children.size() == 0) {
        return true;
    } else {
        
        bool res = true;
        

        if (root->multi_partitions.size() > 4 && root->hybrid_index < root->multi_partitions.size()) {
            
            res = (root->hybrid_index < root->children.size()) || (root->multi_partitions.size() == root->children.size());
            
            std::vector<index_t> hybrid_bucket = root->multi_partitions[root->hybrid_index];
            
            for (index_t hybrid_taxon : hybrid_bucket) {
                candidates.erase(hybrid_taxon);
            }
        }
    
        for (Node * child : root->children) {
            res = res && check_root(child, candidates);
        }
        
        return res;
    }
}

size_t Network::find_child_index_in_children(Node* node, Node * child2find) {
    for (size_t i = 0; i < node->children.size(); i++) {
        if (node->children[i] == child2find) {
            return i;
        }
    }
    return node->children.size();
}

void Network::update_hybird_info(Node* node, index_t old_child_index) {
    
    std::swap(node->multi_partitions[old_child_index], node->multi_partitions.back());
    
    size_t old_root_index = node->multi_partitions.size() - 1;
    
    for (int i = 0; i < 2; i++) {
        if (node->pivots[i] == old_child_index) {
            node->pivots[i] = old_root_index;
        } else if (node->pivots[i] == old_root_index) {
            node->pivots[i] = old_child_index;
        }
    }

    for (auto &kv: node->taxon2partition_id_mapping) {
        if (kv.second == old_child_index) kv.second = old_root_index;
        else if (kv.second == old_root_index) kv.second = old_child_index;
    }

    if (node->hybrid_index == old_child_index) {
        node->hybrid_index = old_root_index;
    } else if (node->hybrid_index == old_root_index) {
        node->hybrid_index = old_child_index;
    }

    for (size_t i = 0; i < node->circle_ordering.size(); i++) {
        if (node->circle_ordering[i] == old_root_index) {
            node->circle_ordering[i] = old_child_index;
        } else if (node->circle_ordering[i] == old_child_index) {
            node->circle_ordering[i] = old_root_index;
        }
    }

}


Node* Network::rerooting_hybrid_tree_above_node(Node * node, Node * child) {
    if (node->parent == NULL) {
        if (node->children.size() >= 2) {
            size_t old_child_index = find_child_index_in_children(node, child);
            update_hybird_info(node, old_child_index);
            std::swap(node->children[old_child_index], node->children.back());
            node->children.pop_back();
            node->parent = child;
            return node;
        } else {
            Node* sib = child->get_sibling();
            sib->parent = child;
            return sib;
        }
    } else {
        Node  *new_child = rerooting_hybrid_tree_above_node(node->parent, node);
        size_t old_root_index = node->multi_partitions.size() - 1;
        size_t old_child_index = find_child_index_in_children(node, child);
        update_hybird_info(node, old_child_index);
        node->children[old_child_index] = new_child;
        node->parent = child;
        return node;
    }
}

Node* Network::reroot_hybrid_tree_on_edge_above_node(Node* node, Tree* tre) {
    Node * new_root = new Node(tre->pseudonym(), false);
    Node * new_child_right = rerooting_hybrid_tree_above_node(node->parent, node);
    new_child_right->parent = new_root;
    node->parent = new_root;
    new_root->children.push_back(node);
    new_root->children.push_back(new_child_right);
    return new_root;

}


std::pair<size_t, std::vector<size_t>> Network::find_above_index_in_cycle(Node *root) {
    std::vector<size_t> circle_ordering_with_sential;
    // not root of the tree
    if (root->multi_partitions.size() == root->children.size() + 1) {
        
        circle_ordering_with_sential.insert(circle_ordering_with_sential.end(), root->circle_ordering.begin(), root->circle_ordering.end());
        
        for (index_t i = 0; i < root->circle_ordering.size(); i++) {
            if (root->circle_ordering[i] == root->children.size()) {
                return {i, circle_ordering_with_sential};
            }
        }
        return {root->circle_ordering.size(), circle_ordering_with_sential};
    }
    // the root of the tree  
    else {
        for (index_t i = 0; i < root->circle_ordering.size() - 1; i++) {
            circle_ordering_with_sential.push_back(root->circle_ordering[i]);
        }
        circle_ordering_with_sential.push_back(-1);
        circle_ordering_with_sential.push_back(root->circle_ordering[root->circle_ordering.size() - 1]);
        return {root->circle_ordering.size() - 1, circle_ordering_with_sential}; // the node is the root of the tree thus no above needed
    }
}

std::string Network::gen_hybrid_id() {
    return "#H" + std::to_string(Network::hybrid_pseudonyms++);
}

NetworkNode *Network::build_network(Node * root, Dict *dict, std::unordered_map<std::string, NetworkNode*> &label2hybrid) {
    if (root->children.size() == 0) {
        return new NetworkNode(dict->index2label(root->index));
    } else {
        // for tree node and unresolved node 
        if (root->multi_partitions.size() <= 4 || root->hybrid_index == root->multi_partitions.size()) {
            // NetworkNode *left = build_network(root->children[0], dict, label2hybrid);
            // NetworkNode *right = build_network(root->children[1], dict, label2hybrid);
            NetworkNode *res = new NetworkNode(dict->index2label(root->index));
            // res->children.push_back(left);
            // res->children.push_back(right);
            // left->parent = res;
            // right->parent = res;
            for (size_t i = 0; i < root->children.size(); i++) {
                NetworkNode *child = build_network(root->children[i], dict, label2hybrid);
                res->children.push_back(child);
                child->parent = res;
            }
            return res;
        }
        // for 4 cycle random  random choose the hybrid 
        // else if (root->multi_partitions.size() == 4) {
            
        // } 
        // for >= 4 cycle resovle it via the circle ordering
        
        else {
            
            auto [above_index_in_cycle_vec, circle_ordering] = find_above_index_in_cycle(root);
            // the current root is not the root of the tree
            
            
            // NetworkNode *current_child_right = new NetworkNode(pseudonym());
            std::string hybrid_name = gen_hybrid_id();
            
            NetworkNode *final_root = new NetworkNode(pseudonym());

            NetworkNode *hybrid_subnet = build_network(root->children[circle_ordering[0]], dict, label2hybrid);
            
            NetworkNode *hybrid_subnet_above = new NetworkNode(hybrid_name);

            label2hybrid[hybrid_name] =  hybrid_subnet_above;

            hybrid_subnet_above->children.push_back(hybrid_subnet);
            
            hybrid_subnet->parent = hybrid_subnet_above;

            NetworkNode *current_left_path = hybrid_subnet_above;
            
           
            
            for (size_t i = 1; i < above_index_in_cycle_vec; i++) {
                
                NetworkNode *left_subnet = build_network(root->children[circle_ordering[i]], dict, label2hybrid);
                
                NetworkNode *next_left_path = new NetworkNode(pseudonym());
                
                next_left_path->children.push_back(left_subnet);
                
                left_subnet->parent = next_left_path;
                
                next_left_path->children.push_back(current_left_path);
                
                current_left_path->parent = next_left_path;
                
                current_left_path = next_left_path;
            }

            final_root->children.push_back(current_left_path);
            current_left_path->parent = final_root;

            NetworkNode *current_right_path = final_root;
           
            
            for (size_t i = above_index_in_cycle_vec + 1; i < circle_ordering.size(); i++) {
                NetworkNode *right_subnet = build_network(root->children[circle_ordering[i]], dict, label2hybrid);
                 
                NetworkNode *next_right_path = new NetworkNode(pseudonym());

                current_right_path->children.push_back(next_right_path);
                next_right_path->parent = current_right_path;

                current_right_path = next_right_path;
                
                current_right_path->children.push_back(right_subnet);
                right_subnet->parent = current_right_path;

            }
            // current_right_path->children.push_back(hybrid_subnet_above);
            // current_right_path->label = hybrid_name;
            // hybrid_subnet_above->hybrid = current_right_path;
            NetworkNode* hybrid_ref = new NetworkNode(hybrid_name);
            hybrid_ref->parent = current_right_path;
            hybrid_ref->hybrid = hybrid_subnet_above; // optional pointer to defining node
            current_right_path->children.push_back(hybrid_ref);
            
            return final_root;
        }
    }
}

NetworkNode *Network::build_network(std::string newick, std::unordered_map<std::string, NetworkNode*> &label2hybrid) {
    if (newick.length() == 0 || newick.at(0) != '(') {
        if (newick[0] == '#') {
            std::size_t sep = newick.find("::");
            std::string lambda = newick.substr(sep + 2, std::string::npos);
            std::string temp = newick.substr(0, sep);
            std::size_t sep_ = temp.find(":");
            std::string length = "";
            if (sep_ != std::string::npos) 
                length = temp.substr(sep_ + 1, std::string::npos);
            std::string label = temp.substr(0, sep_);
            NetworkNode *root = new NetworkNode(label);
            root->lambda = stod__(lambda);
            if (length != "") root->length = stod__(length);
            label2node[label] = root;
            return root;
        }
        else {
            std::size_t sep = newick.find(":");
            std::string length = "";
            if (sep != std::string::npos) 
                length = newick.substr(sep + 1, std::string::npos);
            std::string label = newick.substr(0, sep);
            NetworkNode *root = new NetworkNode(label);
            if (length != "") root->length = stod__(length);
            label2node[label] = root;
            return root;
        }
    }
    else {
        NetworkNode *root = new NetworkNode(pseudonym());
        int k = 1;
        for (int i = 0, j = 0; i < newick.length(); i ++) {
            if (newick.at(i) == '(') j ++;
            if (newick.at(i) == ')') j --;
            if (newick.at(i) == ',' && j == 1) {
                root->children.push_back(build_network(newick.substr(k, i - k), label2hybrid));
                k = i + 1;
            }
        }
        int i = newick.length() - 1;
        while (newick.at(i) != ')') i --;
        root->children.push_back(build_network(newick.substr(k, i - k), label2hybrid));
        std::string branch = newick.substr(i + 1, std::string::npos);
        if (branch.find(";") == std::string::npos) {
            std::size_t sep = branch.find("#");
            if (sep != std::string::npos) {
                //support#hybrid:length::lambda
                //support#hybrid::lambda
                std::string support = branch.substr(0, sep);
                if (support != "") root->support = stod__(support);
                std::size_t sep__ = branch.find("::");
                std::string lambda = branch.substr(sep__ + 2, std::string::npos);
                root->lambda = stod__(lambda);
                std::string temp = branch.substr(sep, sep__);
                std::size_t sep_ = temp.find(":");
                std::string length = "";
                if (sep_ != std::string::npos) 
                    length = temp.substr(sep_ + 1, std::string::npos);
                if (length != "") root->length = stod__(length);
                std::string hybrid = temp.substr(0, sep_);
                label2hybrid[hybrid] = root;
            }
            else {
                //support:length
                sep = branch.find(":");
                std::string support = branch.substr(0, sep);
                if (support != "") root->support = stod__(support);
                std::string length = branch.substr(sep + 1, std::string::npos);
                if (length != "") root->length = stod__(length);
            }
        }
        for (NetworkNode *child : root->children)
            child->parent = root;
        return root;
    }
}

NetworkNode *Network::build_tree_of_blob(NetworkNode *root) {
    if (root->children.size() == 0) {
        NetworkNode *new_root = new NetworkNode(root->label);
        label2node[root->label] = new_root;
        return new_root;
    }
    else {
        NetworkNode *new_root = new NetworkNode(pseudonym());
        for (auto child : root->children) {
            NetworkNode *new_child = build_tree_of_blob(child);
            if (child->bridge) {
                new_root->children.push_back(new_child);
            }
            else {
                for (auto grand_child : new_child->children) 
                    new_root->children.push_back(grand_child);
                new_child->children.clear();
                if (label2node.find(new_child->label) != label2node.end()) 
                    label2node.erase(new_child->label);
                delete new_child;
            }
        }
        for (NetworkNode *child : new_root->children) 
            child->parent = new_root;
        return new_root;
    }
}

NetworkNode *Network::simplify_tree(NetworkNode *root) {
    if (root->children.size() == 0) {
        return root;
    }
    else {
        for (std::size_t i = 0; i < root->children.size(); i ++) 
            root->children[i] = simplify_tree(root->children[i]);
        if (root->children.size() > 1) return root;
        NetworkNode *temp = root->children[0];
        temp->parent = root->parent;
        root->children.clear();
        if (label2node.find(root->label) != label2node.end()) 
            label2node.erase(root->label);
        delete root;
        return temp;
    }
}

std::string Network::display_network(NetworkNode *root) {
    if (root->children.size() == 0) {
        std::string s = root->label + ":" + to_string_p(root->length, 10);
        if (root->label[0] == '#')
            s += "::" + std::to_string(root->lambda);
        return s;
    }
    std::string s = "(";
    for (auto child : root->children) {
        s += display_network(child) + ",";
    }
    s[s.length() - 1] = ')';
    if (root->parent != NULL) {
        // s += std::to_string(root->support);
        if (root->hybrid == NULL) {
            s += ":" + to_string_p(root->length, 10);
        }
        else {
            s += root->hybrid->label;
            s += ":" + to_string_p(root->length, 10);
            s += "::" + to_string_p(root->lambda, 10);
        }
    }
    return s;
}

std::string Network::display_network_basic(NetworkNode *root) {
    if (root->children.size() == 0) 
        return root->label;
    std::string s = "(";
    for (auto child : root->children) 
        s += display_network_basic(child) + ",";
    s[s.length() - 1] = ')';
    // if (root->hybrid != NULL) 
        // s += root->hybrid->label;
    if (!root->label.empty() && root->label.rfind("#H", 0) == 0) {
        s += root->label;
    }
    return s;
}

std::string Network::display_network_bridge(NetworkNode *root) {
    std::string s;
    if (root->children.size() == 0) {
        s = root->label;
    }
    else {
        s = "(";
        for (auto child : root->children) 
            s += display_network_bridge(child) + ",";
        s[s.length() - 1] = ')';
        if (root->hybrid != NULL) 
            s += root->hybrid->label;
    }
    if (root->bridge) s += "^B";
    return s;
}

std::string Network::pseudonym() {
    return "P" + std::to_string(pseudonyms ++);
}

void Network::depth_first_search(NetworkNode *root, NetworkNode *parent, std::size_t depth, std::unordered_map<NetworkNode *, std::pair<std::size_t, std::size_t>> &visited) {
    visited[root] = std::make_pair(depth, depth);
    std::vector<NetworkNode *> descendants;
    if (root->hybrid != NULL) descendants.push_back(root->hybrid);
    if (root->parent != NULL) descendants.push_back(root->parent);
    for (auto child : root->children) descendants.push_back(child);
    for (auto node : descendants) {
        if (node == parent) continue;
        if (visited.find(node) != visited.end()) {
            if (visited[node].first < visited[root].second) 
                visited[root].second = visited[node].first;
        }
        else {
            depth_first_search(node, root, depth + 1, visited);
            if (visited[node].second < visited[root].second)
                visited[root].second = visited[node].second;
        }
    }
    root->bridge = parent != NULL && visited[root].second > visited[parent].first;
}

void Network::compress(double coeffi) {
    compress_length(root, coeffi);
}

void Network::compress_length(NetworkNode *root, double coeffi) {
    root->length *= coeffi;
    for (auto child : root->children) 
        compress_length(child, coeffi);
}

std::vector<NetworkNode *> Network::get_hybridization() {
    return get_hybridization(root);
}

std::vector<NetworkNode *> Network::get_hybridization(NetworkNode *root) {
    std::vector<NetworkNode *> hybrid;
    if (root->hybrid != NULL) {
        hybrid.push_back(root);
    }
    for (auto child : root->children) {
        std::vector<NetworkNode *> sub_hybrid = get_hybridization(child);
        for (auto elem : sub_hybrid) 
            hybrid.push_back(elem);
    }
    return hybrid;
}

std::vector<NetworkNode *> Network::get_blobs() {
    return get_blobs(root);
}

std::vector<NetworkNode *> Network::get_blobs(NetworkNode *root) {
    std::vector<NetworkNode *> blobs;
    if (root->children.size() > 2) 
        blobs.push_back(root);
    for (auto child : root->children) {
        std::vector<NetworkNode *> sub_blobs = get_blobs(child);
        for (auto elem : sub_blobs) 
            blobs.push_back(elem);
    }
    return blobs;
}

std::unordered_map<std::string, NetworkNode*> Network::get_label2hybrid() {
    return label2hybrid;
}
