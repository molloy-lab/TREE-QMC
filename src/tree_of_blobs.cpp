#if ENABLE_TOB
#include "tree.hpp"
#include "rlib_dirs.hpp"
#include "riffle_sort.hpp"

void get_qCFs(std::vector<Tree *> &input, index_t *indices, weight_t *qCFs) {
    index_t temp[4];
    for (index_t i = 0; i < 4; i ++) 
        temp[i] = indices[i];
    std::sort(temp, temp + 4);
    quartet_t q = join(temp);

    weight_t f[3] = {0, 0, 0};
    for (Tree *t : input) {
        index_t topology = t->get_quartet(temp);
        if (topology >= 0) f[topology] += 1;
    }
    qCFs[0] = f[0];
    qCFs[1] = f[1];
    qCFs[2] = f[2];
}




std::pair<Node *, std::vector<index_t>> SpeciesTree::hybrid_info_tree(Node *root, Dict *dict, std::unordered_set<Node *> &false_positive_alpha, std::unordered_set<Node *> &false_positive_beta) {
    std::vector<index_t> minimizers = {};
    if (root->children.size() == 0) {
        // std::cout << "Processing leaf node " << root->index << std::endl;
        Node *new_root = new Node(root->index, false);
        index2node[new_root->index] = new_root;
        return {new_root, minimizers};
    } else {
        
        // std::cout << "Processing node " << root->index << " with " << root->children.size() << " children" << std::endl;
        bool in_false_positive_alpha = false_positive_alpha.find(root) != false_positive_alpha.end();
        
        bool in_false_positive_beta = false_positive_beta.find(root) != false_positive_beta.end();
        
        Node *new_root = new Node(pseudonym(), in_false_positive_alpha || in_false_positive_beta);
        //  Node *new_root = new Node(pseudonym(), in_false_positive_alpha);
        if (in_false_positive_alpha) {
            minimizers.insert(minimizers.end(), std::begin(root->minimizer), std::end(root->minimizer));
        }

        for (Node *child : root->children) {
            // main the child and minimizers that no place to go 
            auto [new_child, child_minimizers] = hybrid_info_tree(child, dict, false_positive_alpha, false_positive_beta);
            
            
            if (!new_child->isfake) {
                new_root->children.push_back(new_child);
                // then the minimizer has a place to go 
                new_child->minimizers = std::move(child_minimizers);
                index2node[new_child->index] = new_child;
            }
            else {
                for (Node *grand_child : new_child->children) 
                    new_root->children.push_back(grand_child);
                
                minimizers.insert(minimizers.end(), child_minimizers.begin(), child_minimizers.end());
                
                new_child->children.clear();
                delete new_child;
            }
        }
        for (Node *new_child : new_root->children) 
            new_child->parent = new_root;
        return {new_root, minimizers};

    }
}

std::vector<std::array<std::vector<index_t>,4 >> sample_quads(const std::vector<std::vector<index_t>>& groups, size_t R, unsigned seed = 42) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> dist(0, 3);
    std::vector<std::array<std::vector<index_t>, 4>> results;
    results.reserve(2*R);

    std::vector<size_t> idx(groups.size());
    std::iota(idx.begin(), idx.end(), 0);

    for (size_t r = 0; r < 2*R; ++r) {
        std::shuffle(idx.begin(), idx.end(), rng);

        std::array<std::vector<index_t>, 4> buckets;
        for (size_t i = 0; i < idx.size(); ++i) {
            int b = static_cast<int>(i % 4);
            const auto& g = groups[idx[i]];
            buckets[b].insert(buckets[b].end(), g.begin(), g.end());
        }
        results.push_back(std::move(buckets));
    }
    return results;
}
// void sample_tree_quads(const std::vector<std::vector<index_t>>& cluster_a, 
//     const std::vector<std::vector<index_t>>& cluster_b, 
//     const std::vector<std::vector<index_t>>&cluster_c,  
//     std::vector<std::array<std::vector<index_t>,4 >>&result,
//     unsigned seed = 42) {

//     std::mt19937 rng(seed);
//     // case 1: root 
//     if (cluster_b.empty() && cluster_c.empty()) {
//         std::vector<size_t> idx(cluster_a.size());
//         std::iota(idx.begin(), idx.end(), 0);
//         std::shuffle(idx.begin(), idx.end(), rng);
//         std::array<std::vector<index_t>, 4> buckets;
//         for (size_t i = 0; i < idx.size(); ++i) {
//             int b = static_cast<int>(i % 4);
//             const auto& g = cluster_a[idx[i]];
//             buckets[b].insert(buckets[b].end(), g.begin(), g.end());
//         }
//         result.push_back(std::move(buckets));
//         std::vector<index_t> above_left = buckets[2];
//         above_left.insert(above_left.end(), buckets[3].begin(), buckets[3].end());
//         std::vector<index_t> above_right = buckets[0];
//         above_right.insert(above_right.end(), buckets[1].begin(), buckets[1].end());
//         sample_tree_quads(buckets[0], buckets[1], above_left);
//         sample_tree_quads(buckets[1], buckets[0], above_right);
//         sample_tree_quads(buckets[2], buckets[3], above_right);
//         sample_tree_quads(buckets[3], buckets[2], above_left);
//     } else if (cluster_a.size() == 1) {
//         return;
//     } else {
//         // shuffle cluster a
//         std::vector<size_t> idx(cluster_a.size());
//         std::iota(idx.begin(), idx.end(), 0);
//         std::shuffle(idx.begin(), idx.end(), rng);
//         std::array<std::vector<index_t>, 2> buckets;
//         for (size_t i = 0; i < idx.size(); ++i) {
//             int b = static_cast<int>(i % 2);
//             const auto& g = cluster_a[idx[i]];
//             buckets[b].insert(buckets[b].end(), g.begin(), g.end());
//         }
//         std::array<std::vector<index_t>, 4> new_quard_partitions;
//         new_quard_partitions[0] = buckets[0];
//         new_quard_partitions[1] = buckets[1];
        
//         for (auto &part : cluster_b) {
//             new_quard_partitions[2].insert(new_quard_partitions[2].end(), part.begin(), part.end());
//         }
        
//         for (auto &part : cluster_c) {
//             new_quard_partitions[3].insert(new_quard_partitions[3].end(), part.begin(), part.end());
//         }
//         // combine cluster b and c
//         result.push_back(std::move(new_quard_partitions));
//     }
//     std::vector<std::array<std::vector<index_t>, 4>> results;
//     results.reserve(2*R);

//     std::vector<size_t> idx(groups.size());
//     std::iota(idx.begin(), idx.end(), 0);

//     for (size_t r = 0; r < 2*R; ++r) {
//         std::shuffle(idx.begin(), idx.end(), rng);

//         std::array<std::vector<index_t>, 4> buckets;
//         for (size_t i = 0; i < idx.size(); ++i) {
//             int b = static_cast<int>(i % 4);
//             const auto& g = groups[idx[i]];
//             buckets[b].insert(buckets[b].end(), g.begin(), g.end());
//         }
//         results.push_back(std::move(buckets));
//     }
//     return results;
// }



// randomly generate minmizers for a blbo node in a tob with expected qcfs
void SpeciesTree::generate_minimizers(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, Node *root, Dict *dict, unsigned long int iter_limit, weight_t alpha) {
    size_t blob_size = root->multi_partitions.size();
    std::vector<std::array<std::vector<index_t>,4 >> quads = sample_quads(root->multi_partitions, blob_size - 3, 42 + root->index);
    size_t count = 0;
    // weight_t current_pvalue = search_quard(gene_trees, quard_vec, current_minimizer);
    
    auto cmp = [](const std::pair<std::array<index_t, 4>, weight_t> &a,
                      const std::pair<std::array<index_t, 4>, weight_t> &b) {
        return a.second > b.second; // min-heap based on weight_t
    };

    std::priority_queue<std::pair<std::array<index_t, 4>, weight_t>,
                        std::vector<std::pair<std::array<index_t, 4>, weight_t>>,
                        decltype(cmp)> pq(cmp);
    for (auto vec : quads) {

        if (count >= blob_size - 3) break;

        index_t current_minimizer[4];
        std::vector<std::vector<index_t>>  quard_vec = {vec[0], vec[1], vec[2], vec[3]};

        weight_t current_pvalue = search_quard_heuristic(qCFs_table, quard_vec, iter_limit, current_minimizer);
        
        if (current_pvalue < alpha) {
            std::cout << "Generated minimizers ";
            for (int i = 0; i < 4; i ++) {
                std::cout << dict->index2label(current_minimizer[i]) << " ";
            }
            std::cout << " with p-value " << current_pvalue << std::endl;

            root->minimizers.push_back(current_minimizer[0]);
            root->minimizers.push_back(current_minimizer[1]);
            root->minimizers.push_back(current_minimizer[2]);
            root->minimizers.push_back(current_minimizer[3]);
            count += 1;
        } else {
            std::cout << "Discarded minimizers ";
            for (int i = 0; i < 4; i ++) {
                std::cout << dict->index2label(current_minimizer[i]) << " ";
            }
            std::cout << " with p-value " << current_pvalue << std::endl;
        }

        pq.push({{current_minimizer[0], current_minimizer[1], current_minimizer[2], current_minimizer[3]}, current_pvalue});
        

    }

    if (count < blob_size - 3) {
        std::cout << "Not enough minimizers found under alpha = " << alpha << ", filling up with best minimizers found." << std::endl;
        while (count < blob_size - 3 && !pq.empty()) {
            auto [minimizer, pvalue] = pq.top();
            pq.pop();

            root->minimizers.push_back(minimizer[0]);
            root->minimizers.push_back(minimizer[1]);
            root->minimizers.push_back(minimizer[2]);
            root->minimizers.push_back(minimizer[3]);
            count += 1;

            std::cout << "Added minimizers ";
            for (int i = 0; i < 4; i ++) {
                std::cout << dict->index2label(minimizer[i]) << " ";
            }
            std::cout << " with p-value " << pvalue << std::endl;
        }
    }
   
}


// randomly generate minmizers for a blbo node in a tob 
void SpeciesTree::generate_minimizers(std::vector<Tree *> &gene_trees, Node *root, Dict *dict, unsigned long int iter_limit, weight_t alpha) {
    size_t blob_size = root->multi_partitions.size();
    std::vector<std::array<std::vector<index_t>,4 >> quads = sample_quads(root->multi_partitions, blob_size - 3, 42 + root->index);
    size_t count = 0;
    // weight_t current_pvalue = search_quard(gene_trees, quard_vec, current_minimizer);
    
    auto cmp = [](const std::pair<std::array<index_t, 4>, weight_t> &a,
                      const std::pair<std::array<index_t, 4>, weight_t> &b) {
        return a.second > b.second; // min-heap based on weight_t
    };

    std::priority_queue<std::pair<std::array<index_t, 4>, weight_t>,
                        std::vector<std::pair<std::array<index_t, 4>, weight_t>>,
                        decltype(cmp)> pq(cmp);
    for (auto vec : quads) {

        if (count >= blob_size - 3) break;

        index_t current_minimizer[4];
        std::vector<std::vector<index_t>>  quard_vec = {vec[0], vec[1], vec[2], vec[3]};
        
        auto get_nodes_from_indices_vec = [&](const std::vector<index_t> &indices) {
            std::vector<Node *> nodes;
            nodes.reserve(indices.size());
            for (index_t idx : indices) {
                nodes.push_back(index2node[idx]);
            }
            return nodes;

        }; 

        std::tuple<std::vector<Node *>, std::vector<Node *>, std::vector<Node *>, std::vector<Node *>> quard_nodes = {
            get_nodes_from_indices_vec(vec[0]),
            get_nodes_from_indices_vec(vec[1]),
            get_nodes_from_indices_vec(vec[2]),
            get_nodes_from_indices_vec(vec[3])
        };
             
        // weight_t current_pvalue = search_quard_heuristic(gene_trees, quard_vec, iter_limit, current_minimizer);
        weight_t current_pvalue = search_3f1a(gene_trees, &quard_nodes, current_minimizer);
        
        auto [_, qcfs] = get_pvalue_and_qCFs(gene_trees, current_minimizer);

        if (current_pvalue < alpha) {
            std::cout << "Generated minimizers ";
            for (int i = 0; i < 4; i ++) {
                std::cout << dict->index2label(current_minimizer[i]) << " ";
            }
            std::cout << " with p-value " << current_pvalue << " with qcfs: ";
            for (const auto& qcf : qcfs) {
                std::cout << qcf << " ";
            }
            std::cout << std::endl;

            root->minimizers.push_back(current_minimizer[0]);
            root->minimizers.push_back(current_minimizer[1]);
            root->minimizers.push_back(current_minimizer[2]);
            root->minimizers.push_back(current_minimizer[3]);
            count += 1;
        } else {
            std::cout << "Discarded minimizers ";
            for (int i = 0; i < 4; i ++) {
                std::cout << dict->index2label(current_minimizer[i]) << " ";
            }
            std::cout << " with p-value " << current_pvalue << " with qcfs: ";
            for (const auto& qcf : qcfs) {
                std::cout << qcf << " ";
            }
            std::cout << std::endl;
        }

        pq.push({{current_minimizer[0], current_minimizer[1], current_minimizer[2], current_minimizer[3]}, current_pvalue});
        

    }

    if (count < blob_size - 3) {
        std::cout << "Not enough minimizers found under alpha = " << alpha << ", filling up with best minimizers found." << std::endl;
        while (count < blob_size - 3 && !pq.empty()) {
            auto [minimizer, pvalue] = pq.top();
            pq.pop();

            root->minimizers.push_back(minimizer[0]);
            root->minimizers.push_back(minimizer[1]);
            root->minimizers.push_back(minimizer[2]);
            root->minimizers.push_back(minimizer[3]);
            count += 1;

            std::cout << "Added minimizers ";
            for (int i = 0; i < 4; i ++) {
                std::cout << dict->index2label(minimizer[i]) << " ";
            }
            std::cout << " with p-value " << pvalue << std::endl;
        }
    }
   
}

// compute taxon to partition mapping for a blob node in a tob using expected qCFS
std::vector<index_t> SpeciesTree::compute_taxon2parition_mapping(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, Node *root, Dict *dict, std::vector<Node *> &hybrid_blob_nodes, std::unordered_set<Node *> & full_leaf_nodes, unsigned long int iter_limit_blob, weight_t alpha) {
    std::vector<index_t> taxon_below;
    if (root->children.size() == 0) {
        taxon_below.push_back(root->index);
    } else {
        index_t partition_id = 0;
        std::unordered_set<index_t > seen_taxa;
        for (Node *child : root->children) {
            std::vector<index_t> child_taxon_below = compute_taxon2parition_mapping(qCFs_table, child, dict, hybrid_blob_nodes, full_leaf_nodes, iter_limit_blob, alpha);
            root->multi_partitions.push_back(child_taxon_below);
            for (index_t taxon : child_taxon_below) {
                taxon_below.push_back(taxon);
                root->taxon2partition_id_mapping[taxon] = partition_id;
            }
            
            seen_taxa.insert(child_taxon_below.begin(), child_taxon_below.end());
            
            partition_id++;
        }

        
        bool flag = true;
        for (Node *taxon : full_leaf_nodes) {
            if (seen_taxa.find(taxon->index) == seen_taxa.end()) {
                if (flag) {
                    root->multi_partitions.push_back({}); // outside taxa partition
                    flag = false;
                }
                root->taxon2partition_id_mapping[taxon->index] = partition_id; // outside taxa
                root->multi_partitions.back().push_back(taxon->index);
            }
        }


        if (root->multi_partitions.size() > 4) {
            hybrid_blob_nodes.push_back(root);

            // for the case that take tob as input so we need to generate minimizers first
            size_t blob_size = root->multi_partitions.size();

            if (root->minimizers.empty()) {
                std::cout << "Generating minimizers for blob id " << root->index << " with " << blob_size << " partitions." << std::endl;
                generate_minimizers(qCFs_table, root, dict, iter_limit_blob, alpha);
            }

            // printing out the minimizers
            std::cout << "Hybrid blob id " << root->index << " with partitions: ";
            for (index_t i = 0; i < root->multi_partitions.size(); i ++) {
                std::cout << "[";
                for (index_t taxon : root->multi_partitions[i]) {
                    std::cout << dict->index2label(taxon) << " ";
                }
                std::cout << "] ";
            }
            std::cout << " and minimizers: ";
            for (index_t taxon : root->minimizers) {
                std::cout << dict->index2label(taxon) << " ";
            }
            std::cout << std::endl;
            std::cout << "the total size of minimizers is " << root->minimizers.size() << std::endl;
        
        }
        
    }

    return taxon_below;
}


std::vector<index_t> SpeciesTree::compute_taxon2parition_mapping(std::vector<Tree *> &gene_trees, Node *root, Dict *dict, std::vector<Node *> &hybrid_blob_nodes, std::unordered_set<Node *> & full_leaf_nodes, unsigned long int iter_limit_blob, weight_t alpha) {
    std::vector<index_t> taxon_below;
    if (root->children.size() == 0) {
        taxon_below.push_back(root->index);
    } else {
        index_t partition_id = 0;
        std::unordered_set<index_t > seen_taxa;
        for (Node *child : root->children) {
            std::vector<index_t> child_taxon_below = compute_taxon2parition_mapping(gene_trees, child, dict, hybrid_blob_nodes, full_leaf_nodes, iter_limit_blob, alpha);
            root->multi_partitions.push_back(child_taxon_below);
            for (index_t taxon : child_taxon_below) {
                taxon_below.push_back(taxon);
                root->taxon2partition_id_mapping[taxon] = partition_id;
            }
            
            seen_taxa.insert(child_taxon_below.begin(), child_taxon_below.end());
            
            partition_id++;
        }

        
        bool flag = true;
        for (Node *taxon : full_leaf_nodes) {
            if (seen_taxa.find(taxon->index) == seen_taxa.end()) {
                if (flag) {
                    root->multi_partitions.push_back({}); // outside taxa partition
                    flag = false;
                }
                root->taxon2partition_id_mapping[taxon->index] = partition_id; // outside taxa
                root->multi_partitions.back().push_back(taxon->index);
            }
        }


        if (root->multi_partitions.size() > 4) {
            hybrid_blob_nodes.push_back(root);

            // for the case that take tob as input so we need to generate minimizers first
            size_t blob_size = root->multi_partitions.size();

            if (root->minimizers.empty()) {
                std::cout << "Generating minimizers for blob id " << root->index << " with " << blob_size << " partitions." << std::endl;
                generate_minimizers(gene_trees, root, dict, iter_limit_blob, alpha);
            }

            // printing out the minimizers
            std::cout << "Hybrid blob id " << root->index << " with partitions: ";
            for (index_t i = 0; i < root->multi_partitions.size(); i ++) {
                std::cout << "[";
                for (index_t taxon : root->multi_partitions[i]) {
                    std::cout << dict->index2label(taxon) << " ";
                }
                std::cout << "] ";
            }
            std::cout << " and minimizers: ";
            for (index_t taxon : root->minimizers) {
                std::cout << dict->index2label(taxon) << " ";
            }
            std::cout << std::endl;
            std::cout << "the total size of minimizers is " << root->minimizers.size() << std::endl;
        
        }
        
    }

    return taxon_below;
}

void SpeciesTree::hybrid_voting(std::vector<Tree *> &gene_trees, Dict *dict, Node* blob_node, unsigned long int iter_limit, std::vector<std::unordered_set<index_t>> &banned_buckets) {

    // add_r_libpaths_and_load(RINS);
    // for (Tree *t : gene_trees) t->LCA_preprocessing();

    std::unordered_map<index_t, int> parititons_votes; 
    parititons_votes.reserve(blob_node->multi_partitions.size());
    std::unordered_set<index_t> banned_buckets_index_of_current_node;
    
    std::vector<std::unordered_set<index_t>> partition_sets;
    partition_sets.reserve(blob_node->multi_partitions.size());

    for (auto &p : blob_node->multi_partitions) {
        partition_sets.emplace_back(p.begin(), p.end());
    }

    for (size_t i = 0; i < banned_buckets.size(); i++) {
        for (size_t j = 0; j < blob_node->multi_partitions.size(); j++) {
            if (std::all_of( banned_buckets[i].begin(),  banned_buckets[i].end(),
            [&](index_t a) {
                return partition_sets[j].count(a) > 0;
            })) {
                banned_buckets_index_of_current_node.insert(j);
            }
        }
    }
    
    for (index_t i = 0; i < blob_node->multi_partitions.size(); i ++) {
        parititons_votes[i] = 0ULL;
        
    }

    for (size_t k = 0; k + 3 < blob_node->minimizers.size(); k += 4) {

        std::unordered_map<index_t, size_t> seen_partitions_in_quad;
        seen_partitions_in_quad.reserve(4);
            
        for (size_t t = 0; t < 4; ++t) {
            index_t minimizer = blob_node->minimizers[k + t];
            
            auto it = blob_node->taxon2partition_id_mapping.find(minimizer);
            index_t partition_id = it->second;

            if (it != blob_node->taxon2partition_id_mapping.end()) {
                seen_partitions_in_quad[partition_id] += 1;
                // seen_partitions_in_quad.insert(partition_id);
            } else {
                seen_partitions_in_quad[partition_id] = 1;
            }
        }

        for (const auto& [partition_id, votes] : seen_partitions_in_quad) {
            if (banned_buckets_index_of_current_node.find(partition_id) != banned_buckets_index_of_current_node.end()) {
                parititons_votes[partition_id] -= 1;
            } else if (votes == 1) {
                parititons_votes[partition_id] += 1;
            } else if (votes >= 2) {
                // parititons_votes[partition_id] = -100;
                parititons_votes[partition_id] += 1;
            }
            
        }
    }

    if (blob_node->multi_partitions.size() < 5) {
        std::cout << "Blob id " << blob_node->index << " has the degree of less than 5, unable to identify the hybridization" << std::endl;
    } else {
            
            // sort the partition votes by values 

            std::vector<std::pair<index_t, int>> sorted_partitions (parititons_votes.begin(), parititons_votes.end());
            std::sort(sorted_partitions.begin(), sorted_partitions.end(), [](auto &a, auto &b) { return a.second > b.second; });
        
            // sorted_partitions.erase(std::remove_if(sorted_partitions.begin(), sorted_partitions.end(), [](std::pair<index_t, int> a) return a.second < 0), sorted_partitions.end());

            // if (sorted_partitions.size() > 5) {
            //     sorted_partitions.resize(5); // keep top 5 partitions only
            // }

            if (sorted_partitions[0].second < 0) {
                std::cout << "Blob id " << blob_node->index << " all partitions conflicted leave it as unresolved. ";
                return;
            }


            
            std::cout << "number of parititons : " << blob_node->multi_partitions.size() << std::endl;
            // print out the top 5 partitions
            std::cout << "Blob id " << blob_node->index << " partitions votes: ";
            for (auto &p : sorted_partitions) {
                std::cout << "Partition " << p.first << " with votes " << p.second << "; ";
            }
            std::cout << std::endl;

            std::cout << "sorted partitions size: " << sorted_partitions.size() << std::endl;


            // iterates all the way of choose 4 from 5 paritions

            if (sorted_partitions[0].second > sorted_partitions[1].second) {
                blob_node->hybrid_index = sorted_partitions[0].first;
            } else {

                weight_t best_pvalue = -1.0;
                index_t best_missing = -1;

                for (index_t missing = 0; missing < 5; ++missing) {
                std::vector<std::vector<index_t>> quard_vec;
                
                
                for (index_t j = 0; j < 5; j++) {
                    if (j != missing) {
                        
                        quard_vec.push_back(blob_node->multi_partitions[sorted_partitions[j].first]);
                        
                    }
                }

                // std::cout << "quard_vec size: " << quard_vec.size() << std::endl;

                
                index_t current_minimizer[4];
                // weight_t current_pvalue = search_quard(gene_trees, quard_vec, current_minimizer);
                weight_t current_pvalue = search_quard_heuristic(gene_trees, quard_vec, iter_limit, current_minimizer);
                // std::cout << "debugging line 210" << std::endl;
                auto [_, qcfs] = get_pvalue_and_qCFs(gene_trees, current_minimizer);
                std::cout << "Trying missing partition " << missing << " with minimizers ";
                for (int i = 0; i < 4; i ++) {
                    std::cout << dict->index2label(current_minimizer[i]) << " ";
                }
                std::cout << " resulting in p-value " << current_pvalue << " with qcfs: ";
                for (const auto& qcf : qcfs) {
                    std::cout << qcf << " ";
                }
                std::cout << std::endl;
                if (best_pvalue < current_pvalue) {
                    best_pvalue = current_pvalue;
                    best_missing = missing;
                }

                
            }
            // std::cout << "best missing : " << best_missing << std::endl;
            // std::cout << "Blob id " << blob_node->index << " identified hybridization partition id: " << sorted_partitions[best_missing].first << " with best p-value: " << best_pvalue << std::endl;
        
            blob_node->hybrid_index = sorted_partitions[best_missing].first;
            
        }

        std::unordered_set<index_t> banned_bucket;
        
        for (size_t j = 0; j < blob_node->multi_partitions.size(); j++) {
            if (j != blob_node->hybrid_index) {
                for (index_t taxon : blob_node->multi_partitions[j]) {
                    banned_bucket.insert(taxon);
                }
            }
        }
        
        banned_buckets.push_back(banned_bucket);
        if (parititons_votes[blob_node->hybrid_index] > 0) {
            
            std::cout << "Hybrid blob id " << blob_node->index << " identified hybridization partition : [";
            
            for (index_t taxon : blob_node->multi_partitions[blob_node->hybrid_index]) {
                std::cout << dict->index2label(taxon) << " ";
            }
            
            std::cout << "]" << std::endl;
        } else {
            blob_node->hybrid_index = blob_node->multi_partitions.size(); // leave it as unresolved
            std::cout << "Blob id " << blob_node->index << " hybridization partition : " << blob_node->hybrid_index << " conflicts with other hybrid node; leave it as unresolved." << std::endl;
        }
        
    }
}

/// hybrid voting with expected qCFS
void SpeciesTree::hybrid_voting(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, Dict *dict, Node* blob_node, unsigned long int iter_limit, std::vector<std::unordered_set<index_t>> &banned_buckets) {

    // add_r_libpaths_and_load(RINS);
    // for (Tree *t : gene_trees) t->LCA_preprocessing();

    std::unordered_map<index_t, int> parititons_votes; 
    parititons_votes.reserve(blob_node->multi_partitions.size());
    std::unordered_set<index_t> banned_buckets_index_of_current_node;
    
    std::vector<std::unordered_set<index_t>> partition_sets;
    partition_sets.reserve(blob_node->multi_partitions.size());

    for (auto &p : blob_node->multi_partitions) {
        partition_sets.emplace_back(p.begin(), p.end());
    }

    for (size_t i = 0; i < banned_buckets.size(); i++) {
        for (size_t j = 0; j < blob_node->multi_partitions.size(); j++) {
            if (std::all_of( banned_buckets[i].begin(),  banned_buckets[i].end(),
            [&](index_t a) {
                return partition_sets[j].count(a) > 0;
            })) {
                banned_buckets_index_of_current_node.insert(j);
            }
        }
    }
    
    for (index_t i = 0; i < blob_node->multi_partitions.size(); i ++) {
        parititons_votes[i] = 0ULL;
        
    }

    for (size_t k = 0; k + 3 < blob_node->minimizers.size(); k += 4) {

        std::unordered_map<index_t, size_t> seen_partitions_in_quad;
        seen_partitions_in_quad.reserve(4);
            
        for (size_t t = 0; t < 4; ++t) {
            index_t minimizer = blob_node->minimizers[k + t];
            
            auto it = blob_node->taxon2partition_id_mapping.find(minimizer);
            index_t partition_id = it->second;

            if (it != blob_node->taxon2partition_id_mapping.end()) {
                seen_partitions_in_quad[partition_id] += 1;
                // seen_partitions_in_quad.insert(partition_id);
            } else {
                seen_partitions_in_quad[partition_id] = 1;
            }
        }

        for (const auto& [partition_id, votes] : seen_partitions_in_quad) {
            if (banned_buckets_index_of_current_node.find(partition_id) != banned_buckets_index_of_current_node.end()) {
                parititons_votes[partition_id] -= 1;
            } else if (votes == 1) {
                parititons_votes[partition_id] += 1;
            } else if (votes >= 2) {
                // parititons_votes[partition_id] = -100;
                parititons_votes[partition_id] += 1;
            }
            
        }
    }

    if (blob_node->multi_partitions.size() < 5) {
        std::cout << "Blob id " << blob_node->index << " has the degree of less than 5, unable to identify the hybridization" << std::endl;
    } else {
            
            // sort the partition votes by values 

        std::vector<std::pair<index_t, int>> sorted_partitions (parititons_votes.begin(), parititons_votes.end());
        std::sort(sorted_partitions.begin(), sorted_partitions.end(), [](auto &a, auto &b) { return a.second > b.second; });
        
        // sorted_partitions.erase(std::remove_if(sorted_partitions.begin(), sorted_partitions.end(), [](std::pair<index_t, int> a) return a.second < 0), sorted_partitions.end());

        // if (sorted_partitions.size() > 5) {
        //     sorted_partitions.resize(5); // keep top 5 partitions only
        // }

        if (sorted_partitions[0].second < 0) {
            std::cout << "Blob id " << blob_node->index << " all partitions conflicted leave it as unresolved. ";
            return;
        }


            
        std::cout << "number of parititons : " << blob_node->multi_partitions.size() << std::endl;
        // print out the top 5 partitions
        std::cout << "Blob id " << blob_node->index << " partitions votes: ";
        for (auto &p : sorted_partitions) {
            std::cout << "Partition " << p.first << " with votes " << p.second << "; ";
        }
        std::cout << std::endl;

        std::cout << "sorted partitions size: " << sorted_partitions.size() << std::endl;


            // iterates all the way of choose 4 from 5 paritions
            
        weight_t best_pvalue = -1.0;
        index_t best_missing = -1;
            
        for (index_t missing = 0; missing < 5; ++missing) {
            std::vector<std::vector<index_t>> quard_vec;
                
                
            for (index_t j = 0; j < 5; j++) {
                if (j != missing) {
                        
                    quard_vec.push_back(blob_node->multi_partitions[sorted_partitions[j].first]);
                        
                }
            }

            // std::cout << "quard_vec size: " << quard_vec.size() << std::endl;

                
            index_t current_minimizer[4];
            // weight_t current_pvalue = search_quard(gene_trees, quard_vec, current_minimizer);
            weight_t current_pvalue = search_quard_heuristic(qCFs_table, quard_vec, iter_limit, current_minimizer);
            // std::cout << "debugging line 210" << std::endl;
            if (best_pvalue < current_pvalue) {
                best_pvalue = current_pvalue;
                best_missing = missing;
            }

                
        }
        // std::cout << "best missing : " << best_missing << std::endl;
        // std::cout << "Blob id " << blob_node->index << " identified hybridization partition id: " << sorted_partitions[best_missing].first << " with best p-value: " << best_pvalue << std::endl;
        
        blob_node->hybrid_index = sorted_partitions[best_missing].first;
        std::unordered_set<index_t> banned_bucket;
        
        for (size_t j = 0; j < blob_node->multi_partitions.size(); j++) {
            if (j != blob_node->hybrid_index) {
                for (index_t taxon : blob_node->multi_partitions[j]) {
                    banned_bucket.insert(taxon);
                }
            }
        }
        
        banned_buckets.push_back(banned_bucket);
        if (parititons_votes[blob_node->hybrid_index] > 0) {
            
            std::cout << "Hybrid blob id " << blob_node->index << " identified hybridization partition : [";
            
            for (index_t taxon : blob_node->multi_partitions[blob_node->hybrid_index]) {
                std::cout << dict->index2label(taxon) << " ";
            }
            
            std::cout << "]" << std::endl;
        } else {
            blob_node->hybrid_index = blob_node->multi_partitions.size(); // leave it as unresolved
            std::cout << "Blob id " << blob_node->index << " hybridization partition : " << blob_node->hybrid_index << " conflicts with other hybrid node; leave it as unresolved." << std::endl;
        }
        
    }
}


void SpeciesTree::pivot_scan(std::vector<Tree *> &gene_trees, Dict *dict, Node *blob_node, unsigned long int iter_limit) {
    
        
    if (blob_node->multi_partitions.size() < 5) {
        std::cout << "Blob id " << blob_node->index << " has the degree of less than 5, unable to identify the hybridization" << std::endl;
        
    } else {
        std::unordered_map<index_t, index_t> pivot_scores;
            
        index_t hybrid_partition_index = blob_node->hybrid_index;
        std::vector<index_t> multi_partitions_index_no_pivot;
            
        for (index_t i = 0; i < blob_node->multi_partitions.size(); i++) {
            pivot_scores[i] = 0;
            if (i != hybrid_partition_index) {
                multi_partitions_index_no_pivot.push_back(i);
            }
        }
        // random shuffle 
        std::random_shuffle(multi_partitions_index_no_pivot.begin(),
                multi_partitions_index_no_pivot.end());

        // index_t partition_dmmy = multi_partitions_index_no_pivot[0];
        
        // for (index_t fixed = 1;  fixed < multi_partitions_index_no_pivot.size(); fixed ++) {
            
        //     for (index_t k = fixed + 1; k < multi_partitions_index_no_pivot.size(); k ++) {
                
        //     index_t partition_k = multi_partitions_index_no_pivot[k];
        //     index_t partititon_fixed = multi_partitions_index_no_pivot[fixed];
        //     // below is the minimum p-value method
        //     // std::vector<std::vector<index_t>> quad = {blob_node->multi_partitions[partititon_fixed], 
        //     //     blob_node->multi_partitions[partition_dmmy],
        //     //     blob_node->multi_partitions[partition_k],
        //     //     blob_node->multi_partitions[hybrid_partition_index]
        //     // };

        //     // index_t minimizer[4];
        //     // weight_t min_pvalue = search_quard_heuristic(gene_trees, quad, iter_limit, minimizer);
        //     // index_t hybrid_taxon_index_pos = -1;
        //     // for (index_t i = 0; i < 4; i++) {
        //     //     if (blob_node->taxon2partition_id_mapping[minimizer[i]] == hybrid_partition_index) {
        //     //         hybrid_taxon_index_pos = i;
        //     //         break;
        //     //     }
        //     // }
        //     // index_t hybrid_taxon_index = minimizer[hybrid_taxon_index_pos];
        //     // std::array<std::array<index_t,4>,2> displayed_quartets = computed_displayed_quartet_toplogy(minimizer);
        //     // std::array<index_t, 2> hybrid_siblings = siblings_in_two_best_topologies(displayed_quartets, hybrid_taxon_index);
        //     //above is the minimum p-value method

        //     // comment out below is average qCFs method
        //     std::array<index_t, 4> quad_ids = {
        //         multi_partitions_index_no_pivot[partititon_fixed],
        //         partition_dmmy,
        //         partition_k,
        //         hybrid_partition_index
        //     };

        //     std::sort(quad_ids.begin(), quad_ids.end());
            

        //     std::array<weight_t, 3> qCFs = freq_three_toplogies(gene_trees, blob_node, quad_ids, dict->size());
        //     auto sib_parts = hybrid_siblings_from_top2_qcfs(qCFs, quad_ids, hybrid_partition_index);
        //     std::array<index_t, 2> hybrid_siblings = {sib_parts[0], sib_parts[1]};
        //     // comment out above is average qCFs method
            


        //     if (hybrid_siblings[0] == partititon_fixed || hybrid_siblings[1] == partititon_fixed) {
        //         pivot_scores[partititon_fixed] += 1;
        //     }
        //     }
        
        // }
        // // get the top 2 highest pivot scores
        // index_t highest_pivot = -1;
        // index_t second_highest_pivot = -1;
        // index_t highest_score = -1;
        // index_t second_highest_score = -1;
        // for (const auto& [partition_id, score] : pivot_scores) {
        //     if (score > highest_score) {
        //         second_highest_score = highest_score;
        //         second_highest_pivot = highest_pivot;
        //         highest_score = score;
        //         highest_pivot = partition_id;
        //     } else if (score > second_highest_score) {
        //         second_highest_score = score;
        //         second_highest_pivot = partition_id;
        //     }
        // }
        // blob_node->pivots[0] = highest_pivot;
        // blob_node->pivots[1] = second_highest_pivot;

        // commened out the pivot scan below 
        index_t partititon_i = multi_partitions_index_no_pivot[0];
        index_t partittion_j = multi_partitions_index_no_pivot[1];
        for (index_t k = 2; k < multi_partitions_index_no_pivot.size(); k ++) {
                
            index_t partition_k = multi_partitions_index_no_pivot[k];

            std::array<index_t, 4> quad_ids = {
                partititon_i,
                partittion_j,
                partition_k,
                hybrid_partition_index
            };
                

            std::array<weight_t, 3> qCFs = freq_three_toplogies(gene_trees, blob_node, quad_ids, dict->size());
            auto sib_parts = hybrid_siblings_from_top2_qcfs(qCFs, quad_ids, hybrid_partition_index);
            std::array<index_t, 2> hybrid_siblings = {sib_parts[0], sib_parts[1]};
            partititon_i = hybrid_siblings[0];
            partittion_j = hybrid_siblings[1];

        }

        
        blob_node->pivots[0] = partititon_i < partittion_j ? partititon_i : partittion_j;
        blob_node->pivots[1] = partititon_i < partittion_j ? partittion_j : partititon_i;
        // commened out the pivot scan above

        std::cout << "blob id : " << blob_node->index << " has hybrid bucket: [ ";
        for (index_t taxon : blob_node->multi_partitions[blob_node->hybrid_index]) {
            std::cout << dict->index2label(taxon) << " ";
        }
        std::cout << "] and pivot bucket 0 : [ ";
            
        for (index_t taxon : blob_node->multi_partitions[blob_node->pivots[0]]) {
            std::cout << dict->index2label(taxon) << " ";
        }
        std::cout << "] " << std::endl;

        std::cout << "] and pivot bucket 1 : [ ";
            
        for (index_t taxon : blob_node->multi_partitions[blob_node->pivots[1]]) {
            std::cout << dict->index2label(taxon) << " ";
        }
        std::cout << "] " << std::endl;


        // Print all pivot scores for all buckets
        // std::cout << "pivot scores by bucket id: ";
        // for (index_t pid = 0; pid < (index_t)blob_node->multi_partitions.size(); ++pid) {
            // std::cout << pid << "=" << pivot_scores[pid];
            // if (pid + 1 < (index_t)blob_node->multi_partitions.size()) std::cout << ", ";
        // }
        std::cout << std::endl;
        std::cout << "piviot candidate 1 = " << blob_node->pivots[0] << std::endl;
        std::cout << "piviot candidate 2 = " << blob_node->pivots[1] << std::endl;

    }

}

std::array<index_t, 2> SpeciesTree::hybrid_siblings_from_top2_qcfs(const std::array<weight_t, 3>& qCFs,
                               const std::array<index_t,4>& quad_ids, index_t hybrid_partition_id) {
    int hpos = -1;
    for (int i = 0; i < 4; ++i) {
        if (quad_ids[i] == hybrid_partition_id) { hpos = i; break; }
    }
    if (hpos < 0) {
        std::cout << "Error: hybrid_partition_id not found in quad_ids." << std::endl;
        return { (index_t)-1, (index_t)-1 };
    }

    // Helper: for a given topology t in {0,1,2}, return the sibling POSITION of hpos.
    auto sibling_pos = [&](int topo) -> int {
        // topo 0: (0,1)|(2,3)
        if (topo == 0) {
            if (hpos == 0) return 1;
            if (hpos == 1) return 0;
            if (hpos == 2) return 3;
            return 2; // hpos == 3
        }
        // topo 1: (0,2)|(1,3)
        if (topo == 1) {
            if (hpos == 0) return 2;
            if (hpos == 2) return 0;
            if (hpos == 1) return 3;
            return 1; // hpos == 3
        }
        // topo 2: (0,3)|(1,2)
        // pairs: (0,3) and (1,2)
        if (hpos == 0) return 3;
        if (hpos == 3) return 0;
        if (hpos == 1) return 2;
        return 1; // hpos == 2
    };

    // pick top-2 topologies by qCF (same tie-break)
    int best = 0, second = 1;
    auto better = [&](int a, int b) -> bool {
        if (qCFs[a] > qCFs[b]) return true;
        if (qCFs[a] < qCFs[b]) return false;
        return a < b;
    };

    if (better(1, best)) { second = best; best = 1; }
    else { second = 1; }

    if (better(2, best)) {
        second = best;
        best = 2;
    } else if (better(2, second)) {
        second = 2;
    }

    // sibling PARTITION IDS (indices into multi_partitions)
    int sibpos1 = sibling_pos(best);
    int sibpos2 = sibling_pos(second);

    return { quad_ids[sibpos1], quad_ids[sibpos2] };
}


bool SpeciesTree::is_bucket_i_less_than_bucket_j(index_t partition_i, index_t partition_j, index_t pivot_index, Node* blob_node, std::vector<Tree *> gene_trees, unsigned long int iter_limit, size_t &failed_counts) {
    index_t hybrid_index = blob_node->hybrid_index;

    // iter_limit = std::ceil(2 * iter_limit / std::log(dict->size()));
    // below is minmimum p-value method
    // std::vector<std::vector<index_t>> quad = {blob_node->multi_partitions[partition_i], 
    //                 blob_node->multi_partitions[partition_j],
    //                 blob_node->multi_partitions[pivot_index],
    //                 blob_node->multi_partitions[hybrid_index]
    //             };
    
    // index_t minimizer[4];

    // weight_t min_pvalue = search_quard_heuristic(gene_trees, quad, iter_limit, minimizer);
    
    
    // index_t hybrid_taxon_index_pos = -1;
    
    // for (index_t i = 0; i < 4; i++) {
    //     if (blob_node->taxon2partition_id_mapping[minimizer[i]] == hybrid_index) {
    //         hybrid_taxon_index_pos = i;
    //         break;
    //     }
    // }
    
    // index_t hybrid_taxon_index = minimizer[hybrid_taxon_index_pos];

    // std::array<std::array<index_t,4>,2> displayed_quartets = computed_displayed_quartet_toplogy(minimizer);

    // std::array<index_t, 2> hybrid_siblings = siblings_in_two_best_topologies(displayed_quartets, hybrid_taxon_index);
    // above is minmimum p-value method

    // comment out below is average qCFs method
    std::array<index_t, 4> quad_ids = {partition_i, partition_j, pivot_index, hybrid_index};
    // std::cout << "current 4 tuple: ";
    // std::cout << " [ "
    // for (index_t i : blob_node->multi_partitions[partition_i]) {
    //         std::cout << dict->index2label(i) << " ";
    //     }
    // std::cout << " ] | [ ";
    // for (index_t i : blob_node->multi_partitions[partition_j]) {
    //         std::cout << dict->index2label(i) << " ";
    //     }
    // std::cout << " ] | [ ";
    // for (index_t i : blob_node->multi_partitions[pivot_index]) {

    //         std::cout << dict->index2label(i) << " ";
    //     }
    // std::cout << " ] | [ ";
    // for (index_t i : blob_node->multi_partitions[hybrid_index]) {
    //         std::cout << dict->index2label(i) << " ";
    //     }
    // std::cout << " ] ";

    // std::sort(quad_ids.begin(), quad_ids.end());
    std::array<weight_t, 3> qCFs = freq_three_toplogies(gene_trees, blob_node, quad_ids, dict->size());
    // std::cout << "with qCFs: ";
    // for (index_t t = 0; t < 3; t++) {
    //     std::cout << qCFs[t] << " ";
    // }
    // std::cout << std::endl;

    auto sib_parts = hybrid_siblings_from_top2_qcfs(qCFs, quad_ids, hybrid_index);
    std::array<index_t, 2> hybrid_siblings = {sib_parts[0], sib_parts[1]};
    // comment out above is average qCFs method

    bool seen_pivot_with_hybrid = false;
    bool bucket_i_less_than_j = true;
    for (index_t i = 0; i < 2; i++) {
        if (hybrid_siblings[i] == pivot_index) {
            seen_pivot_with_hybrid = true;
        }
    }

    if (!seen_pivot_with_hybrid) {
        // std::cout << "Warnning for blob id : " << blob_node->index << " pivot : " << pivot_index <<" bucket : " << partition_i << " and bucket " << partition_j << " We did not see hybrid|pivot siblings in all displayed quartet topology" << std::endl;
        failed_counts += 1;

        std::cout << "failed 4-tuple and hybrid sibling and pvalue: ";
        for (index_t i : blob_node->multi_partitions[partition_i]) {
            std::cout << dict->index2label(i) << " ";
        }
        std::cout << " | ";
        for (index_t i : blob_node->multi_partitions[partition_j]) {
            std::cout << dict->index2label(i) << " ";
        }
        std::cout << " | ";
        for (index_t i : blob_node->multi_partitions[pivot_index]) {
            std::cout << dict->index2label(i) << " ";
        }
        std::cout << " | ";
        for (index_t i : blob_node->multi_partitions[hybrid_index]) {
            std::cout << dict->index2label(i) << " ";
        }; 
        
        std::cout << " | hybrid: [";
        for (auto taxon : blob_node->multi_partitions[hybrid_index]) {
            std::cout << dict->index2label(taxon) << " ";
        }
        std::cout << "] ";
        std::cout << " hybrid siblings: ";
        for (index_t t = 0; t < 2; t++) {
            std::cout << "[ ";
            for (auto taxon : blob_node->multi_partitions[hybrid_siblings[t]]) {
                std::cout << dict->index2label(taxon) << " ";
            }
            std::cout << "] ";
        }
        std::cout << " | qCFs: ";
        for (index_t t = 0; t < 3; t++) {
            std::cout << qCFs[t] << " ";
        }
        std::cout << std::endl;
    }

    for (index_t i = 0; i < 2; i++) {
        if (hybrid_siblings[i] == pivot_index) {
            continue;
        } else if (hybrid_siblings[i] == partition_i) {
            bucket_i_less_than_j = false;
        }
    }
    // std::cout << "finished comparing bucket " << partition_i << " and bucket " << partition_j << " with pivot " << pivot_index << ". Result: " << (bucket_i_less_than_j ? "true" : "false") << std::endl;
    return bucket_i_less_than_j;


}

// using expected qCFs
bool SpeciesTree::is_bucket_i_less_than_bucket_j(index_t partition_i, index_t partition_j, index_t pivot_index, Node* blob_node, std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, unsigned long int iter_limit, size_t &failed_counts) {
    index_t hybrid_index = blob_node->hybrid_index;

    // iter_limit = std::ceil(2 * iter_limit / std::log(dict->size()));
    // below is minmimum p-value method
    std::vector<std::vector<index_t>> quad = {blob_node->multi_partitions[partition_i], 
                    blob_node->multi_partitions[partition_j],
                    blob_node->multi_partitions[pivot_index],
                    blob_node->multi_partitions[hybrid_index]
                };
    

    index_t minimizer[4] = {quad[0][0], quad[1][0], quad[2][0], quad[3][0]};

    std::sort(minimizer, minimizer + 4);
    // std::cout << "quad taxon: ";
    // for (index_t t = 0; t < 4; t++) {
    //     std::cout << dict->index2label(minimizer[t]) << " ";
    // }
    auto [min_pvalue, qcf] = get_pvalue_and_qCFs(qCFs_table, minimizer);
    // std::cout << "min pvalue: " << min_pvalue << std::endl;
    // std::cout << "qCFs: ";
    // for (index_t t = 0; t < 3; t++) {
    //     std::cout << qcf[t] << " ";
    // }
    // std::cout << std::endl;

    index_t hybrid_taxon_index_pos = -1;
    
    for (index_t i = 0; i < 4; i++) {
        if (blob_node->taxon2partition_id_mapping[minimizer[i]] == hybrid_index) {
            hybrid_taxon_index_pos = i;
            break;
        }
    }
    
    index_t hybrid_taxon_index = minimizer[hybrid_taxon_index_pos];

    std::array<std::array<index_t,4>,2> displayed_quartets = computed_displayed_quartet_toplogy(minimizer, qcf);

    std::array<index_t, 2> hybrid_siblings = siblings_in_two_best_topologies(displayed_quartets, hybrid_taxon_index);
    // above is minmimum p-value method
    // std::cout << "hybrid siblings: ";
    
    // for (index_t t = 0; t < 2; t++) {
    //     std::cout << dict->index2label(hybrid_siblings[t]) << " ";
    // }
    // std::cout << std::endl;

    bool seen_pivot_with_hybrid = false;
    bool bucket_i_less_than_j = true;
    for (index_t i = 0; i < 2; i++) {
        if (blob_node->taxon2partition_id_mapping[hybrid_siblings[i]] == pivot_index) {
            seen_pivot_with_hybrid = true;
        }
    }

    if (!seen_pivot_with_hybrid) {
        // std::cout << "Warnning for blob id : " << blob_node->index << " pivot : " << pivot_index <<" bucket : " << partition_i << " and bucket " << partition_j << " We did not see hybrid|pivot siblings in all displayed quartet topology" << std::endl;
        failed_counts += 1;
        std::cout << "failed 4-tuple and hybrid sibling and pvalue: ";
        for (index_t t = 0; t < 4; t++) {
            std::cout << dict->index2label(minimizer[t]) << " ";
        }
        std::cout << " | pvalue: " << min_pvalue;
        std::cout << " | hybrid: [";
        for (auto taxon : blob_node->multi_partitions[hybrid_index]) {
            std::cout << dict->index2label(taxon) << " ";
        }
        std::cout << "] ";
        std::cout << " hybrid siblings: ";
        for (index_t t = 0; t < 2; t++) {
            std::cout << dict->index2label(hybrid_siblings[t]) << " ";
        }
        std::cout << " | qCFs: ";
        for (index_t t = 0; t < 3; t++) {
            std::cout << qcf[t] << " ";
        }
        std::cout << std::endl;

    }

    for (index_t i = 0; i < 2; i++) {
        if (blob_node->taxon2partition_id_mapping[hybrid_siblings[i]] == pivot_index) {
            continue;
        } else if (blob_node->taxon2partition_id_mapping[hybrid_siblings[i]] == partition_i) {
            bucket_i_less_than_j = false;
        }
    }
    // std::cout << "finished comparing bucket " << partition_i << " and bucket " << partition_j << " with pivot " << pivot_index << ". Result: " << (bucket_i_less_than_j ? "true" : "false") << std::endl;
    return bucket_i_less_than_j;


}


void SpeciesTree::circle_sorting_enmuerate_pivots(std::vector<Tree *> &gene_trees, unsigned long int iter_limit, Node * blob_node) {

    if (blob_node->multi_partitions.size() == 4) {
        std::cout << "blob id : " << blob_node->index << " has the degree of less than 5, unable to identify the hybridization " << std::endl;
        return; 
    } 

    const index_t hybrid_partition_index = blob_node->hybrid_index;

    // Build base list of candidate partitions (excluding hybrid and pivot later)
    std::vector<index_t> base_ids;
    base_ids.reserve(blob_node->multi_partitions.size());
    for (index_t pid = 0; pid < blob_node->multi_partitions.size(); ++pid) {
        if (pid == hybrid_partition_index) continue;
        base_ids.push_back(pid);
    }

        

    // Track best attempt
    size_t best_failed = std::numeric_limits<size_t>::max();
    index_t best_pivot = blob_node->multi_partitions.size();
    std::vector<index_t> best_circle;

    auto try_with_pivot = [&](index_t pivot_partition_index) -> void {

        // ids excluding pivot & hybrid
        std::vector<index_t> ids;
        ids.reserve(blob_node->multi_partitions.size());
        for (index_t pid : base_ids) {
            if (pid == pivot_partition_index) continue;
            ids.push_back(pid);
        }

        size_t failed_counts = 0;

        std::sort(ids.begin(), ids.end(),
                    [&](index_t a, index_t b) {
                        return is_bucket_i_less_than_bucket_j(a, b, pivot_partition_index, blob_node, 
                                                            gene_trees, iter_limit,
                                                            failed_counts);
                    });

        // auto comp = [&](index_t a, index_t b) {
        //                 return is_bucket_i_less_than_bucket_j(a, b, pivot_partition_index, blob_node, 
        //                                                     gene_trees, iter_limit,
        //                                                     failed_counts);
        //             };


        // riffleSort(ids, 0.0, comp);

        std::vector<index_t> circle_ordering;
        circle_ordering.reserve(blob_node->multi_partitions.size());
        circle_ordering.push_back(hybrid_partition_index);
        circle_ordering.push_back(pivot_partition_index);
        circle_ordering.insert(circle_ordering.end(), ids.begin(), ids.end());

        // Keep best (min failed). Tie-break: keep the first found.
        if (best_pivot == blob_node->multi_partitions.size() || failed_counts < best_failed) {
            best_failed = failed_counts;
            best_pivot = pivot_partition_index;
            best_circle = std::move(circle_ordering);
        }
    };

    // 1) Try the first pivot

    for (auto pivot_partition_index : base_ids) {
        std::cout << "Trying pivot partition id : " << pivot_partition_index << std::endl;
        try_with_pivot(pivot_partition_index);
        std::cout << "Finished trying pivot partition id : " << pivot_partition_index << " with failed counts : " << best_failed << std::endl;
        if (best_failed == 0) break; // early stop if perfect

    }
    

        // Apply best result
    if (best_pivot != blob_node->multi_partitions.size()) {
        blob_node->circle_ordering = best_circle;

        std::cout << "blob id " << blob_node->index
                    << " chosen pivot=" << best_pivot
                    << " failed_count=" << best_failed
                    << " circle ordering:\n";

        for (index_t pid : blob_node->circle_ordering) {
            std::cout << "  bucket " << pid << ": [ ";
            if (pid >= 0 && pid < (index_t)blob_node->multi_partitions.size()) {
                for (index_t taxon : blob_node->multi_partitions[pid]) {
                    std::cout << dict->index2label(taxon) << " ";
                }
            } else {
                std::cout << "(invalid bucket id) ";
            }
            std::cout << "]\n";
        }
    } else {
        std::cout << "Warning: blob id " << blob_node->index
                    << " no valid pivot found for circle sorting.\n";
    }
}

// circle sorting with enumerating pivots using expected qCFS
void SpeciesTree::circle_sorting_enmuerate_pivots(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, unsigned long int iter_limit, Node * blob_node) {

    if (blob_node->multi_partitions.size() == 4) {
        std::cout << "blob id : " << blob_node->index << " has the degree of less than 5, unable to identify the hybridization " << std::endl;
        return; 
    } 

    const index_t hybrid_partition_index = blob_node->hybrid_index;

    // Build base list of candidate partitions (excluding hybrid and pivot later)
    std::vector<index_t> base_ids;
    base_ids.reserve(blob_node->multi_partitions.size());
    for (index_t pid = 0; pid < blob_node->multi_partitions.size(); ++pid) {
        if (pid == hybrid_partition_index) continue;
        base_ids.push_back(pid);
    }

        

    // Track best attempt
    size_t best_failed = std::numeric_limits<size_t>::max();
    index_t best_pivot = blob_node->multi_partitions.size();
    std::vector<index_t> best_circle;

    auto try_with_pivot = [&](index_t pivot_partition_index) -> void {

        // ids excluding pivot & hybrid
        std::vector<index_t> ids;
        ids.reserve(blob_node->multi_partitions.size());
        for (index_t pid : base_ids) {
            if (pid == pivot_partition_index) continue;
            ids.push_back(pid);
        }

        size_t failed_counts = 0;

        std::sort(ids.begin(), ids.end(),
                    [&](index_t a, index_t b) {
                        return is_bucket_i_less_than_bucket_j(a, b, pivot_partition_index, blob_node, 
                                                            qCFs_table, iter_limit,
                                                            failed_counts);
                    });

        std::vector<index_t> circle_ordering;
        circle_ordering.reserve(blob_node->multi_partitions.size());
        circle_ordering.push_back(hybrid_partition_index);
        circle_ordering.push_back(pivot_partition_index);
        circle_ordering.insert(circle_ordering.end(), ids.begin(), ids.end());

        // Keep best (min failed). Tie-break: keep the first found.
        if (best_pivot == blob_node->multi_partitions.size() || failed_counts <= best_failed) {
            best_failed = failed_counts;
            best_pivot = pivot_partition_index;
            best_circle = std::move(circle_ordering);
        }
    };

    // 1) Try the first pivot

    for (auto pivot_partition_index : base_ids) {
        std::cout << "Trying pivot partition id : " << pivot_partition_index << std::endl;
        try_with_pivot(pivot_partition_index);
        std::cout << "Finished trying pivot partition id : " << pivot_partition_index << " with failed counts : " << best_failed << std::endl;
        if (best_failed == 0) break; // early stop if perfect

    }
    

        // Apply best result
    if (best_pivot != blob_node->multi_partitions.size()) {
        blob_node->circle_ordering = best_circle;

        std::cout << "blob id " << blob_node->index
                    << " chosen pivot=" << best_pivot
                    << " failed_count=" << best_failed
                    << " circle ordering:\n";

        for (index_t pid : blob_node->circle_ordering) {
            std::cout << "  bucket " << pid << ": [ ";
            if (pid >= 0 && pid < (index_t)blob_node->multi_partitions.size()) {
                for (index_t taxon : blob_node->multi_partitions[pid]) {
                    std::cout << dict->index2label(taxon) << " ";
                }
            } else {
                std::cout << "(invalid bucket id) ";
            }
            std::cout << "]\n";
        }
    } else {
        std::cout << "Warning: blob id " << blob_node->index
                    << " no valid pivot found for circle sorting.\n";
    }
}

void SpeciesTree::circle_sorting(std::vector<Tree *> &gene_trees, unsigned long int iter_limit, Node * blob_node) {

    if (blob_node->multi_partitions.size() == 4) {
        std::cout << "blob id : " << blob_node->index << " has the degree of less than 5, unable to identify the hybridization " << std::endl;
        // std::cout << "Using default circle ordering for 4 partitions. " << std::endl;
        // blob_node->circle_ordering.push_back(0);
        // blob_node->circle_ordering.push_back(1);
        // blob_node->circle_ordering.push_back(3);
        // blob_node->circle_ordering.push_back(2);
            
        return; 
    } 

    const index_t hybrid_partition_index = blob_node->hybrid_index;

    // Build base list of candidate partitions (excluding hybrid and pivot later)
    std::vector<index_t> base_ids;
    base_ids.reserve(blob_node->multi_partitions.size());
    for (index_t pid = 0; pid < blob_node->multi_partitions.size(); ++pid) {
        if (pid == hybrid_partition_index) continue;
        base_ids.push_back(pid);
    }

        

    // Track best attempt
    size_t best_failed = std::numeric_limits<size_t>::max();
    index_t best_pivot = -1;
    std::vector<index_t> best_circle;

    auto try_with_pivot = [&](index_t pivot_partition_index) -> void {

        // ids excluding pivot & hybrid
        std::vector<index_t> ids;
        ids.reserve(blob_node->multi_partitions.size());
        for (index_t pid : base_ids) {
            if (pid == pivot_partition_index) continue;
            ids.push_back(pid);
        }

        size_t failed_counts = 0;

        std::sort(ids.begin(), ids.end(),
                    [&](index_t a, index_t b) {
                        return is_bucket_i_less_than_bucket_j(a, b, pivot_partition_index, blob_node, 
                                                            gene_trees, iter_limit,
                                                            failed_counts);
                    });

        // auto comp = [&](index_t a, index_t b) {
        //                 return is_bucket_i_less_than_bucket_j(a, b, pivot_partition_index, blob_node, 
        //                                                     gene_trees, iter_limit,
        //                                                     failed_counts);
        //             };


        // riffleSort(ids, 0.0, comp);

        std::vector<index_t> circle_ordering;
        circle_ordering.reserve(blob_node->multi_partitions.size());
        circle_ordering.push_back(hybrid_partition_index);
        circle_ordering.push_back(pivot_partition_index);
        circle_ordering.insert(circle_ordering.end(), ids.begin(), ids.end());

        // Keep best (min failed). Tie-break: keep the first found.
        if (failed_counts < best_failed) {
            best_failed = failed_counts;
            best_pivot = pivot_partition_index;
            best_circle = std::move(circle_ordering);
        }
    };

    // 1) Try the first pivot
    try_with_pivot(blob_node->pivots[0]);

    // 2) If not perfect, try other pivots and pick min failed
    if (best_failed > 0) {
        for (size_t pi = 1; pi < 2; ++pi) {
            try_with_pivot(blob_node->pivots[pi]);
            if (best_failed == 0) break; // early stop if perfect
        }
    }

        // Apply best result
    if (best_pivot >= 0) {
        blob_node->circle_ordering = best_circle;

        std::cout << "blob id " << blob_node->index
                    << " chosen pivot=" << best_pivot
                    << " failed_count=" << best_failed
                    << " circle ordering:\n";

        for (index_t pid : blob_node->circle_ordering) {
            std::cout << "  bucket " << pid << ": [ ";
            if (pid >= 0 && pid < (index_t)blob_node->multi_partitions.size()) {
                for (index_t taxon : blob_node->multi_partitions[pid]) {
                    std::cout << dict->index2label(taxon) << " ";
                }
            } else {
                std::cout << "(invalid bucket id) ";
            }
            std::cout << "]\n";
        }
    } else {
        std::cout << "Warning: blob id " << blob_node->index
                    << " no valid pivot found for circle sorting.\n";
    }
}
  

// void SpeciesTree::circle_sorting(std::vector<Tree *> gene_trees, size_t iter_limit, std::vector<Node *> &hybrid_blob_nodes) {

//     for (Node *blob_node : hybrid_blob_nodes) {
//         if (blob_node->multi_partitions.size() > 4) {
//             std::vector<index_t> multi_partitions_index_no_pivot_hybrid;
//             index_t hybrid_partition_index = blob_node->hybrid_index;
//             index_t pivot_partition_index = blob_node->pivots[0];
//             std::vector<index_t> multi_partitions_index_no_pivot;
//             std::vector<index_t> circle_ordering;
//             circle_ordering.push_back(hybrid_partition_index);
//             circle_ordering.push_back(pivot_partition_index);
//             for (index_t i = 0; i < blob_node->multi_partitions.size(); i++) {
                
//                 if (i != hybrid_partition_index && i != pivot_partition_index) {
//                     multi_partitions_index_no_pivot.push_back(i);
//                 }
//             }

//             size_t failed_counts = 0;
//             std::sort(multi_partitions_index_no_pivot.begin(), multi_partitions_index_no_pivot.end(), 
//                 [&](index_t a, index_t b) {
//                 return is_bucket_i_less_than_bucket_j(a, b, blob_node, gene_trees, iter_limit, failed_counts);
//             });

//             circle_ordering.insert(circle_ordering.end(), multi_partitions_index_no_pivot.begin(), multi_partitions_index_no_pivot.end());

//             blob_node->circle_ordering = circle_ordering;

//             if (failed_counts > 0) {
//                 size_t failed_counts_another = 0;

//             }
            
//             std::cout << "blob id " << blob_node->index << " circle ordering:\n";
//             for (index_t pid : blob_node->circle_ordering) {
//                 std::cout << "  bucket " << pid << ": [ ";
//                 if (pid >= 0 && pid < (index_t)blob_node->multi_partitions.size()) {
//                     for (index_t taxon : blob_node->multi_partitions[pid]) {
//                         std::cout << dict->index2label(taxon) << " ";
//                     }
//                 } else {
//                     std::cout << "(invalid bucket id) ";
//                 }
//                 std::cout << "]\n";
//             }
//         }
        
//     }

// }


// compute network with expected qcfs table
SpeciesTree::SpeciesTree(Tree *input, Dict *dict, weight_t alpha, weight_t beta, std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, unsigned long int iter_limit_blob) {
        
        std::cout << "Contracting branches with alpha = " << alpha << " and beta = " << beta << std::endl;
        
        this->dict = dict;
        
        std::vector<Node *> internal;
        
        std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> bips;
        
        input->get_bipartitions(&internal, &bips);
        
        std::cout << bips.size() << " branches to test" << std::endl;
        
        std::unordered_set<Node *> false_positive_beta;
        
        std::unordered_set<Node *> false_positive_alpha;

        // std::unordered_set<Node *> false_positive;

        for (index_t i = 0; i < internal.size(); i ++) {
                
            // if (internal[i]->min_pvalue < alpha) {
            //     false_positive_alpha.insert(internal[i]);
            //     // false_positive.insert(internal[i]);
            // } else if (internal[i]->max_pvalue > beta) {
            //     false_positive_beta.insert(internal[i]);
            //     // false_positive.insert(internal[i]);
            // }
                
        }

        std::cout << "number of branches in the input tree: " << internal.size() << std::endl;

        auto [new_root, root_min] = hybrid_info_tree(input->root, dict, false_positive_alpha, false_positive_beta);
        root = new_root;
        root->minimizers = root_min;
        index2node[root->index] = root;
        std::cout << "Number of branches to contract under alpha criterion: " << false_positive_alpha.size() << std::endl;
        std::cout << "Number of branches to contract under beta criterion: " << false_positive_beta.size() << std::endl;
        

        std::cout << "The tree of blobs : " << input->display_tree_basic(root) << std::endl;
        std::unordered_set<Node *> full_leaf_nodes;
        get_leaf_set(root, &full_leaf_nodes);
        std::unordered_set<Node *> full_leaf_indices;
        for (Node *leaf : full_leaf_nodes) {
            full_leaf_indices.insert(leaf);
        }


        add_r_libpaths_and_load(RINS);

        std::vector<Node *> hybrid_blob_nodes;

        compute_taxon2parition_mapping(qCFs_table, root, dict, hybrid_blob_nodes, full_leaf_indices, iter_limit_blob, alpha);

        // for debuging hybrid_infor_tree
        
        // root = build_refinement(input->root, false_positive);
        // std::cout << "TOB using build_refinement: " << input->display_tree_basic(root) << std::endl;

        

        std::sort(hybrid_blob_nodes.begin(), hybrid_blob_nodes.end(), [](const Node * a, Node * b){
            return a->multi_partitions.size() > b->multi_partitions.size();
        });
        std::vector<std::unordered_set<index_t>> banned_buckets;

        

        std::cout << "Printing output number of nontrival blobs:" << hybrid_blob_nodes.size() << std::endl;
        
        for (Node * blob_node : hybrid_blob_nodes) {
            qCFs_average_cache.clear();
            hybrid_voting(qCFs_table, dict, blob_node, iter_limit_blob, banned_buckets);
            std::cout << "Done with hybrid voting for blob id : " << blob_node->index << std::endl;
            
            
            if (blob_node->hybrid_index == blob_node->multi_partitions.size()) {
                std::cout << "blob id : " << blob_node->index << " has no hybridization detected, skip circle sorting " << std::endl;
                continue;
            }

            // pivot_scan(gene_trees, dict, blob_node, iter_limit_blob);
            // circle_sorting(gene_trees, iter_limit_blob, blob_node);

            circle_sorting_enmuerate_pivots(qCFs_table, iter_limit_blob, blob_node);
            std::cout << "Done with circle sorting for blob id : " << blob_node->index << std::endl;

        }
        
        
    }

// given tob-annotatiation tree and compute the m-partion maps and its inverse for all  nontrival blobs and the hybrization parititons
SpeciesTree::SpeciesTree(Tree *input, Dict *dict, weight_t alpha, weight_t beta, std::vector<Tree *> &gene_trees, unsigned long int iter_limit_blob) {
        
        std::cout << "Contracting branches with alpha = " << alpha << " and beta = " << beta << std::endl;
        
        this->dict = dict;
        
        std::vector<Node *> internal;
        
        std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> bips;
        
        input->get_bipartitions(&internal, &bips);
        
        std::cout << bips.size() << " branches to test" << std::endl;
        
        std::unordered_set<Node *> false_positive_beta;
        
        std::unordered_set<Node *> false_positive_alpha;

        // std::unordered_set<Node *> false_positive;

        for (index_t i = 0; i < internal.size(); i ++) {
                
            if (internal[i]->min_pvalue < alpha) {
                false_positive_alpha.insert(internal[i]);
                // false_positive.insert(internal[i]);
            } 
            // else if (internal[i]->max_pvalue > beta) {
                // false_positive_beta.insert(internal[i]);
                // false_positive.insert(internal[i]);
            // }
                
        }

        // std::unordered_set<Node*> internal_set(internal.begin(), internal.end());

        // std::function<void(Node*)> show_missing = [&](Node* u){
        //     if (!u) return;
        //     if (!u->children.empty() && !internal_set.count(u)) {
        //     std::cerr << "[MISSING] ptr=" << u
        //             << " deg=" << u->children.size()
        //             << " min_p=" << u->min_pvalue
        //             << " max_p=" << u->max_pvalue
        //             << "\n";
        //         }
        //     for (auto* c : u->children) show_missing(c);
        //     };
        // show_missing(input->root);

        std::cout << "number of branches in the input tree: " << internal.size() << std::endl;

        auto [new_root, root_min] = hybrid_info_tree(input->root, dict, false_positive_alpha, false_positive_beta);
        root = new_root;
        root->minimizers = root_min;
        index2node[root->index] = root;
        std::cout << "Number of branches to contract under alpha criterion: " << false_positive_alpha.size() << std::endl;
        std::cout << "Number of branches to contract under beta criterion: " << false_positive_beta.size() << std::endl;
        

        std::cout << "The tree of blobs : " << input->display_tree_basic(root) << std::endl;
        std::unordered_set<Node *> full_leaf_nodes;
        get_leaf_set(root, &full_leaf_nodes);
        std::unordered_set<Node *> full_leaf_indices;
        for (Node *leaf : full_leaf_nodes) {
            full_leaf_indices.insert(leaf);
        }


        add_r_libpaths_and_load(RINS);
        for (Tree *t : gene_trees) t->LCA_preprocessing();

        std::vector<Node *> hybrid_blob_nodes;

        compute_taxon2parition_mapping(gene_trees, root, dict, hybrid_blob_nodes, full_leaf_indices, iter_limit_blob, alpha);

        // for debuging hybrid_infor_tree
        
        // root = build_refinement(input->root, false_positive);
        // std::cout << "TOB using build_refinement: " << input->display_tree_basic(root) << std::endl;

        

        std::sort(hybrid_blob_nodes.begin(), hybrid_blob_nodes.end(), [](const Node * a, Node * b){
            return a->multi_partitions.size() > b->multi_partitions.size();
        });
        std::vector<std::unordered_set<index_t>> banned_buckets;

        

        std::cout << "Printing output number of nontrival blobs:" << hybrid_blob_nodes.size() << std::endl;
        
        for (Node * blob_node : hybrid_blob_nodes) {
            qCFs_average_cache.clear();
            hybrid_voting(gene_trees, dict, blob_node, iter_limit_blob, banned_buckets);
            std::cout << "Done with hybrid voting for blob id : " << blob_node->index << std::endl;
            
            
            if (blob_node->hybrid_index == blob_node->multi_partitions.size()) {
                std::cout << "blob id : " << blob_node->index << " has no hybridization detected, skip circle sorting " << std::endl;
                continue;
            }

            // pivot_scan(gene_trees, dict, blob_node, iter_limit_blob);
            // circle_sorting(gene_trees, iter_limit_blob, blob_node);

            circle_sorting_enmuerate_pivots(gene_trees, iter_limit_blob, blob_node);
            std::cout << "Done with circle sorting for blob id : " << blob_node->index << std::endl;

        }
        
        
    }


// load p-value and just testing it get tob
SpeciesTree::SpeciesTree(Tree *input, Dict *dict, weight_t alpha, weight_t beta, bool enable_split_test) {
    std::cout << "Contracting branches with alpha = " << alpha << " and beta = " << beta << std::endl;
    if (enable_split_test) {
        std::cout << "Also consider branches with nonzero split match and mismatch count as false positive for refinement" << std::endl;
    }
    this->dict = dict;
    
    std::vector<Node *> internal;
    std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> bips;
    
    input->get_bipartitions(&internal, &bips);
    
    std::cout << bips.size() << " branches to test" << std::endl;
    std::unordered_set<Node *> false_positive;
    
    for (index_t i = 0; i < internal.size(); i ++) {
        if (internal[i]->min_pvalue < alpha || internal[i]->max_pvalue > beta) {
            false_positive.insert(internal[i]);
        } else if (enable_split_test && internal[i]->split_match_count > 0 && internal[i]->split_mismatch_count > 0) {
            std::cout << "Branch id " << i << " has nonzero split match and mismatch count, consider it as false positive for refinement. split match count: " << internal[i]->split_match_count << " split mismatch count: " << internal[i]->split_mismatch_count << std::endl;
            false_positive.insert(internal[i]);
        }
    }
    
    root = build_refinement(input->root, false_positive);
}
///////////above is for net-cs //////////////////////////////////


// 2f2a search algorithm O(n^3)
SpeciesTree::SpeciesTree(std::vector<Tree *> &input, Dict *dict, SpeciesTree* display) {
    std::cout << "Constructing tree of blobs using 2-fix-2-alter search" << std::endl;

    add_r_libpaths_and_load(RINS);
    for (Tree *t : input) t->LCA_preprocessing();
    this->dict = display->dict;
    display->refine();

    //std::string mode = "n";
    //display->annotate(input, mode);

    std::vector<Node *> internal;
    std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> bips;
    display->get_bipartitions(&internal, &bips);
    std::cout << bips.size() << " branches to test" << std::endl;
    std::unordered_set<Node *> false_positive;

    for (std::size_t i = 0; i < internal.size(); i ++) {
        std::cout << "Testing branch id " << i << ", ";
        internal[i]->blob_id = i;

        if (internal[i]->isfake) {
            false_positive.insert(internal[i]);
            std::cout << "fake ***" << std::endl;
            continue;
        }

        // Search for min p-value for quartet tree tests
        weight_t min;
        index_t minimizer[4];
        size_t split_match_count = 0;
        size_t split_mismatch_count = 0;

        std:: cout << "testing branch with split : [";
        for (Node *n : bips[i].first) {
            std::cout << dict->index2label(n->index) << ",";
        }
        std::cout << "]" << std::endl << "and [";
        for (Node *n : bips[i].second) {
            std::cout << dict->index2label(n->index) << ",";
        }
        std::cout << "]" << std::endl;


        min = search_2f2a(input, bips[i].first, bips[i].second, minimizer, split_match_count, split_mismatch_count);
        internal[i]->min_pvalue = min;

        // Get qCFs that yielded the min p-value
        weight_t min_f[3];
        get_qCFs(input, minimizer, min_f);
        internal[i]->min_f[0] = min_f[0];
        internal[i]->min_f[1] = min_f[1];
        internal[i]->min_f[2] = min_f[2];
        
        internal[i]->minimizer[0] = minimizer[0];
        internal[i]->minimizer[1] = minimizer[1];
        internal[i]->minimizer[2] = minimizer[2];
        internal[i]->minimizer[3] = minimizer[3];
        internal[i]->split_match_count = split_match_count;
        internal[i]->split_mismatch_count = split_mismatch_count;

        // Appy quartet star test
        weight_t max = -1.0;
        if ((min_f[0] + min_f[1] + min_f[2]) > 0)
            max = pvalue_star(min_f);
        //if ((internal[i]->f[0] + internal[i]->f[1] + internal[i]->f[2]) > 0)
        //    max = pvalue_star(internal[i]->f);
        //if (iter_limit != 0)
        //    max = search_star(input, bips[i].first, bips[i].second, iter_limit);
        //else
        //    max = search_star(input, bips[i].first, bips[i].second);
        internal[i]->max_pvalue = max;

        // Write to standard out
        std::cout << "QTT: " << min << "; ";
        std::cout << "QST: " << max << "; ";
        std::cout << "qCF: [" << min_f[0] << "/" << min_f[1] << "/" << min_f[2] << "]; ";
        std::cout << "minimizer: [" << dict->index2label(minimizer[0]) << "/" << dict->index2label(minimizer[1]) << "/" << dict->index2label(minimizer[2]) << "/" << dict->index2label(minimizer[3]) << "] ";
        std::cout << "split match/mismatch count: [" << split_match_count << "/" << split_mismatch_count << "]";
        std::cout << std::endl;

        //std::cout << "QTT: " << min << " ";
        //std::cout << "qCF: [" << internal[i]->min_f[0] << "/" << internal[i]->min_f[1] << "/" << internal[i]->min_f[2] << "] ";
        //std::cout << "minimizer: [" << dict->index2label(minimizer[0]) << "/" << dict->index2label(minimizer[1]) << "/" << dict->index2label(minimizer[2]) << "/" << dict->index2label(minimizer[3]) << "] ";
        //std::cout << "QST: " << max << " ";
        //std::cout << "[" << internal[i]->f[0] << "/" << internal[i]->f[1] << "/" << internal[i]->f[2] << "]";
        //std::cout << std::endl;

        // Write to a table...
    }
    if (display->root->children.size() == 2) 
        false_positive.insert(display->root->children[1]);

    root = build_refinement(display->root, false_positive);
}

// 3 fix 1 alter search algorithm O(n^2)
SpeciesTree::SpeciesTree(std::vector<Tree *> &input, Dict *dict, 
                         SpeciesTree* display, 
                         unsigned long int iter_limit_blob, 
                         bool three_fix_one_alter,
                         bool two_fix_two_alter, 
                         bool is_quard) {
    
    if (!three_fix_one_alter && !is_quard && !two_fix_two_alter) {
       SpeciesTree(input, dict, display, iter_limit_blob);
       return;
    }

    if (two_fix_two_alter && three_fix_one_alter) {
        std::cerr << "Error: both two_fix_two_alter and three_fix_one_alter cannot be true at the same time." << std::endl;
        return;
    }

    if (two_fix_two_alter && is_quard) {
        std::cerr << "Error: both two_fix_two_alter and is_quard cannot be true at the same time." << std::endl;
        return;
    }

    if (two_fix_two_alter) {
        SpeciesTree(input, dict, display);
        return;
    }
    
    if (is_quard){
        std::cout << "Constructing tree of blobs using quard search" << std::endl;
    } else if (three_fix_one_alter) {
        std::cout << "Constructing tree of blobs using 3-fix-1-alter search" << std::endl;
    }

    add_r_libpaths_and_load(RINS);
    for (Tree *t : input) t->LCA_preprocessing();

    display->refine();

    this->dict = display->dict;

    //std::string mode = "n";
    //display->annotate(input, mode);
    
    std::vector<Node *> internal;
    
    std::vector<std::tuple<std::vector<Node *>, std::vector<Node *>, std::vector<Node *>, std::vector<Node *>>> quads;
    display->get_quardpartitions(&internal, &quads, dict);
    std::cout << quads.size() << " branches to test" << std::endl;

    std::unordered_set<Node *> false_positive;
    for (index_t i = 0; i < internal.size(); i ++) {
        std::cout << "Testing branch id " << i << ", ";
        internal[i]->blob_id = i;

        if (internal[i]->isfake) {
            false_positive.insert(internal[i]);
            std::cout << "fake ***" << std::endl;
            continue;
        }


        // Search for min p-value for quartet tree tests
        weight_t min;
        index_t minimizer[4];
        size_t match_count = 0;
        size_t mismatch_count = 0;
        if (three_fix_one_alter) {
            
            min = search_3f1a(input, &quads[i], minimizer);
        }   
        else if (is_quard)
            min = search_quard(input, &quads[i], minimizer);
        internal[i]->min_pvalue = min;

        // Get qCFs that yielded the min p-value
        weight_t min_f[3];
        get_qCFs(input, minimizer, min_f);
        internal[i]->min_f[0] = min_f[0];
        internal[i]->min_f[1] = min_f[1];
        internal[i]->min_f[2] = min_f[2];

        internal[i]->minimizer[0] = minimizer[0];
        internal[i]->minimizer[1] = minimizer[1];
        internal[i]->minimizer[2] = minimizer[2];
        internal[i]->minimizer[3] = minimizer[3];
        if (two_fix_two_alter) {
            internal[i]->split_match_count = match_count;
            internal[i]->split_mismatch_count = mismatch_count;
        }
        // Apply quartet star test
        weight_t max = -1.0;
        if ((min_f[0] + min_f[1] + min_f[2]) > 0)
            max = pvalue_star(min_f);
        internal[i]->max_pvalue = max;
        //if ((internal[i]->f[0] + internal[i]->f[1] + internal[i]->f[2]) > 0)
        //    max = pvalue_star(internal[i]->f);

        // Write to standard out
        std::cout << "QTT: " << min << "; ";
        std::cout << "QST: " << max << "; ";
        std::cout << "qCF: [" << min_f[0] << "/" << min_f[1] << "/" << min_f[2] << "]; ";
        std::cout << "minimizer: [" << dict->index2label(minimizer[0]) << "/" << dict->index2label(minimizer[1]) << "/" << dict->index2label(minimizer[2]) << "/" << dict->index2label(minimizer[3]) << "]";
        std::cout << std::endl;
        //std::cout << "QTT: " << min << " ";
        //std::cout << "qCF: [" << internal[i]->min_f[0] << "/" << internal[i]->min_f[1] << "/" << internal[i]->min_f[2] << "] ";
        //std::cout << "minimizer: [" << dict->index2label(minimizer[0]) << "/" << dict->index2label(minimizer[1]) << "/" << dict->index2label(minimizer[2]) << "/" << dict->index2label(minimizer[3]) << "] ";
        //std::cout << "QST: " << max << " ";
        //std::cout << "[" << internal[i]->f[0] << "/" << internal[i]->f[1] << "/" << internal[i]->f[2] << "]";
        //std::cout << std::endl;
    }

    if (display->root->children.size() == 2) {
        false_positive.insert(display->root->children[1]);
    }

    this->dict = display->dict;
    root = build_refinement(display->root, false_positive);
}


// O(n*cn*klogn), O(kn^3logn) if c is O(n)
SpeciesTree::SpeciesTree(std::vector<Tree *> &input, Dict *dict, 
                         SpeciesTree* display,
                         unsigned long int iter_limit_blob) {
    std::cout << "Constructing tree of blobs" << std::endl;
    add_r_libpaths_and_load(RINS);
    for (Tree *t : input) t->LCA_preprocessing();
    this->dict = display->dict;
    display->refine();

    //std::string mode = "n";
    //display->annotate(input, mode);

    std::vector<Node *> internal;
    std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> bips;
    display->get_bipartitions(&internal, &bips);
    std::cout << bips.size() << " branches to test" << std::endl;
    std::unordered_set<Node *> false_positive;
    size_t iter_limit = iter_limit_blob;

    for (std::size_t i = 0; i < internal.size(); i ++) {
        std::cout << "Testing branch id " << i << ", ";
        internal[i]->blob_id = i;

        if (internal[i]->isfake) {
            false_positive.insert(internal[i]);
            std::cout << "fake ***" << std::endl;
            continue;
        }

        // Search for min p-value for quartet tree tests
        weight_t min;
        index_t minimizer[4];
        if (iter_limit != 0)
            min = search(input, bips[i].first, bips[i].second, iter_limit, minimizer);
        else
            min = search(input, bips[i].first, bips[i].second, minimizer);
        internal[i]->min_pvalue = min;

        // Get qCFs that yielded the min p-value
        weight_t min_f[3];
        get_qCFs(input, minimizer, min_f);
        internal[i]->min_f[0] = min_f[0];
        internal[i]->min_f[1] = min_f[1];
        internal[i]->min_f[2] = min_f[2];
        
        internal[i]->minimizer[0] = minimizer[0];
        internal[i]->minimizer[1] = minimizer[1];
        internal[i]->minimizer[2] = minimizer[2];
        internal[i]->minimizer[3] = minimizer[3];
    

        // Appy quartet star test
        weight_t max = -1.0;
        if ((min_f[0] + min_f[1] + min_f[2]) > 0)
            max = pvalue_star(min_f);
        //if ((internal[i]->f[0] + internal[i]->f[1] + internal[i]->f[2]) > 0)
        //    max = pvalue_star(internal[i]->f);
        //if (iter_limit != 0)
        //    max = search_star(input, bips[i].first, bips[i].second, iter_limit);
        //else
        //    max = search_star(input, bips[i].first, bips[i].second);
        internal[i]->max_pvalue = max;

        // Write to standard out
        std::cout << "QTT: " << min << "; ";
        std::cout << "QST: " << max << "; ";
        std::cout << "qCF: [" << min_f[0] << "/" << min_f[1] << "/" << min_f[2] << "]; ";
        std::cout << "minimizer: [" << dict->index2label(minimizer[0]) << "/" << dict->index2label(minimizer[1]) << "/" << dict->index2label(minimizer[2]) << "/" << dict->index2label(minimizer[3]) << "] ";
        std::cout << std::endl;

        //std::cout << "QTT: " << min << " ";
        //std::cout << "qCF: [" << internal[i]->min_f[0] << "/" << internal[i]->min_f[1] << "/" << internal[i]->min_f[2] << "] ";
        //std::cout << "minimizer: [" << dict->index2label(minimizer[0]) << "/" << dict->index2label(minimizer[1]) << "/" << dict->index2label(minimizer[2]) << "/" << dict->index2label(minimizer[3]) << "] ";
        //std::cout << "QST: " << max << " ";
        //std::cout << "[" << internal[i]->f[0] << "/" << internal[i]->f[1] << "/" << internal[i]->f[2] << "]";
        //std::cout << std::endl;

        // Write to a table...
    }
    if (display->root->children.size() == 2) 
        false_positive.insert(display->root->children[1]);

    root = build_refinement(display->root, false_positive);
}


std::string SpeciesTree::to_string_pvalue() {
    return display_tree_pvalue(root) + ";";
}

std::string SpeciesTree::display_tree_pvalue(Node *root) {
    if (root->children.size() == 0) 
        return dict->index2label(root->index);
    std::string s = "(";
    for (Node * node : root->children) 
        s += display_tree_pvalue(node) + ",";
    s[s.size() - 1] = ')';
    if (root->parent != NULL && (root->parent->parent != NULL || (root->parent->parent == NULL && root == root->parent->children[0]))) {
        std::ostringstream ss;
        // serialize the minimizer and qCFs
        ss << std::scientific << std::setprecision(12) 
                              << "'[" 
                              << "blob_id=" << std::to_string(root->blob_id)
                              << ";qtt_p=" << (double) root->min_pvalue
                              << ";qst_p=" << (double) root->max_pvalue 
                              << ";qcf_1=" << std::to_string((int) root->min_f[0])
                              << ";qcf_2=" << std::to_string((int) root->min_f[1])
                              << ";qcf_3=" << std::to_string((int) root->min_f[2])
                              << ";minimizer=" << dict->index2label(root->minimizer[0]) << "/" << dict->index2label(root->minimizer[1]) << "/" << dict->index2label(root->minimizer[2]) << "/" << dict->index2label(root->minimizer[3])
                             << ";split_match_count=" << root->split_match_count
                             << ";split_mismatch_count=" << root->split_mismatch_count
                              << "]'";
        return s + ss.str();
    }
    else {
        return s;
    }
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

void Tree::get_quardpartition(Node *root, std::vector<Node *> *A, std::vector<Node *> *B, std::vector<Node *> *C, std::vector<Node *> *D, Dict *dict) {
    //std::cout << "Getting quardpartition for node " << root->index << std::endl;

    std::vector<Node *> leaves;
    get_leaves(this->root, &leaves);

    if (root->children.size() < 2) {
            std::cerr << "Error: less than two children" << std::endl;
            exit(1);
        }

    if (root->children.size() > 2) {
        std::cerr << "Error: more than two children" << std::endl;
        exit(1);
    }
    
    Node *c_1, *c_2, *s_1;

    c_1 = root->children[0];
    c_2 = root->children[1];

    if (root != this->root && root->parent != this->root) {
        s_1 = root->parent->children[0] == root ? root->parent->children[1] : root->parent->children[0];

    } else if (root != this->root && root->parent == this->root){
        s_1 = (this->root->children[0] == root) ? this->root->children[1]->children[0] : this->root->children[0]->children[0];
    }

    std::unordered_set<Node *> leaf_set1;
    std::unordered_set<Node *> leaf_set2;
    std::unordered_set<Node *> leaf_set3;
    get_leaf_set(c_1, &leaf_set1);
    get_leaf_set(c_2, &leaf_set2);
    get_leaf_set(s_1, &leaf_set3);
    for (Node *leaf : leaves) {
        if (leaf_set1.find(leaf) != leaf_set1.end()) 
            A->push_back(leaf);
        else if (leaf_set2.find(leaf) != leaf_set2.end()) 
            B->push_back(leaf);
        else if (leaf_set3.find(leaf) != leaf_set3.end()) 
            C->push_back(leaf);
        else 
            D->push_back(leaf);
    }
    // std::string quard_str = display_quardpartitions(*A, *B, *C, *D, dict);
    // // checking sum of sizes 
    // if (A->size() + B->size() + C->size() + D->size() != leaves.size()) {
    //     std::cerr << "Error in quardpartition: sizes do not add up!" << std::endl;
    //     exit(1);
    // }
    // std::cout << "Quardpartition: " << quard_str << std::endl;
    
}

void Tree::get_quardpartitions(Node *root, std::vector<Node *> *internal, std::vector<std::tuple<std::vector<Node *>, std::vector<Node *>, std::vector<Node *>, std::vector<Node *>>> *quards, std::unordered_set<uint64_t>* seen_edges, Dict *dict) {
    
    if (root->children.size() == 0) return ;

    
    if (root->parent != NULL && (root->parent->parent != NULL || root->parent->children.size() > 2 || root == root->parent->children[0])) {
        std::vector<Node *> A, B, C, D;
        get_quardpartition(root, &A, &B, &C, &D, dict);
        auto quard = std::make_tuple(A, B, C, D); 
        internal->push_back(root); // std::cout << "?" << root->support << std::endl;
        quards->push_back(quard);
    }
    
    for (Node *child : root->children) {
        get_quardpartitions(child, internal, quards, seen_edges, dict);
    }
}

void Tree::get_quardpartitions(std::vector<Node *> *internal, 
                               std::vector<std::tuple<std::vector<Node *>,
                               std::vector<Node *>, std::vector<Node *>,std::vector<Node *>>> *quards, Dict *dict) {
    std::unordered_set<uint64_t> seen_edges;
    get_quardpartitions(root, internal, quards, &seen_edges, dict);
}


std::string Tree::display_quardpartitions(std::vector<Node *> &A, std::vector<Node *> &B, std::vector<Node *> &C, std::vector<Node *> &D, Dict *dict) {
    std::string s = "";
    for (Node *a : A) s += dict->index2label(a->index) + ",";
    s = s.substr(0, s.size() - 1);
    s += "| ";
    for (Node *b : B) s += dict->index2label(b->index) + ",";
    s = s.substr(0, s.size() - 1);
    s += "| ";
    for (Node *c : C) s += dict->index2label(c->index) + ",";
    s = s.substr(0, s.size() - 1);
    s += "| ";
    for (Node *d : D) s += dict->index2label(d->index) + ",";
    s = s.substr(0, s.size() - 1);
    return s;
}


weight_t SpeciesTree::search(std::vector<Tree *> &input, 
                             std::vector<Node *> &A, 
                             std::vector<Node *> &B, 
                             index_t* minimizer) {
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

                    weight_t score = get_pvalue(input, temp);
                    if (min < 0 || score < min) {
                        min = score;
                        minimizer[0] = temp[0]; 
                        minimizer[1] = temp[1]; 
                        minimizer[2] = temp[2]; 
                        minimizer[3] = temp[3];
                    }
                }
            }
        }
    }
    //std::cout << "naive iter: " << count << std::endl;
    return min;
}


weight_t SpeciesTree::search_quard(std::vector<Tree *> &input, std::tuple<std::vector<Node *>, std::vector<Node *>, std::vector<Node *>, std::vector<Node *>> *quad, index_t* minimizer) {
    index_t i[4] = {0, 0, 0, 0};
    weight_t min = -1;
    auto& t = *quad;

    for (i[0] = 0; i[0] < std::get<0>(t).size(); i[0] ++) {
        for (i[1] = 0; i[1] < std::get<1>(t).size(); i[1] ++) {
            for (i[2] = 0; i[2] < std::get<2>(t).size(); i[2] ++) {
                for (i[3] = 0; i[3] < std::get<3>(t).size(); i[3] ++) {
                    index_t temp[4];
                    temp[0] = std::get<0>(t)[i[0]]->index;
                    temp[1] = std::get<1>(t)[i[1]]->index;
                    temp[2] = std::get<2>(t)[i[2]]->index;
                    temp[3] = std::get<3>(t)[i[3]]->index;
                    weight_t score = get_pvalue(input, temp);
                    if (min < 0 || score < min) {
                        min = score;
                        minimizer[0] = temp[0]; minimizer[1] = temp[1]; minimizer[2] = temp[2]; minimizer[3] = temp[3];
                    }
                }
            }
        }
    }

    return min;
}


weight_t SpeciesTree::search_quard(std::vector<Tree *> &input, std::vector<std::vector<index_t>> &quad, index_t* minimizer) {
    index_t i[4] = {0, 0, 0, 0};
    weight_t min = -1;

    for (i[0] = 0; i[0] < quad[0].size(); i[0] ++) {
        for (i[1] = 0; i[1] < quad[1].size(); i[1] ++) {
            for (i[2] = 0; i[2] < quad[2].size(); i[2] ++) {
                for (i[3] = 0; i[3] < quad[3].size(); i[3] ++) {
                    index_t temp[4];
                    temp[0] = quad[0][i[0]];
                    temp[1] = quad[1][i[1]];
                    temp[2] = quad[2][i[2]];
                    temp[3] = quad[3][i[3]];
                    // std::cout << "debugging line 755 " << std::endl;
                    weight_t score = get_pvalue(input, temp);
                    // std::cout << "done with get_pvalue " << std::endl;
                    if (min < 0 || score < min) {
                        min = score;
                        minimizer[0] = temp[0]; minimizer[1] = temp[1]; minimizer[2] = temp[2]; minimizer[3] = temp[3];
                    }
                }
            }
        }
    }

    return min;
}

/// for average qcfs

void SpeciesTree::postorder_nodes(Node* root, std::vector<Node*>& out) {
    out.clear();
    if (!root) return;

    std::vector<Node*> st1;
    std::vector<Node*> st2;

    st1.push_back(root);
    while (!st1.empty()) {
        Node* u = st1.back();
        st1.pop_back();
        st2.push_back(u);
        for (Node* c : u->children) {
            if (c) st1.push_back(c);
        }
    }
    while (!st2.empty()) {
        out.push_back(st2.back());
        st2.pop_back();
    }
}

static inline double I_of_C(const std::array<std::array<std::size_t,4>,3>& C) {
    const auto& A = C[0];
    const auto& B = C[1];
    const auto& R = C[2];
    return
        (double)A[0]*B[1]*R[2]*R[3] + (double)A[1]*B[0]*R[2]*R[3] +
        (double)A[2]*B[3]*R[0]*R[1] + (double)A[3]*B[2]*R[0]*R[1] +
        (double)R[0]*B[1]*A[2]*A[3] + (double)R[1]*B[0]*A[2]*A[3] +
        (double)R[2]*B[3]*A[0]*A[1] + (double)R[3]*B[2]*A[0]*A[1] +
        (double)A[0]*R[1]*B[2]*B[3] + (double)A[1]*R[0]*B[2]*B[3] +
        (double)A[2]*R[3]*B[0]*B[1] + (double)A[3]*R[2]*B[0]*B[1];
}

// overload: buckets is std::array<vec,4>
static inline void build_bucket_of(std::vector<int8_t>& bucket_of,
                                   const std::array<std::vector<index_t>,4>& buckets,
                                   index_t num_taxa)
{
    bucket_of.assign((size_t)num_taxa, (int8_t)-1);
    for (int b = 0; b < 4; ++b) {
        for (index_t t : buckets[b]) {
            if ((size_t)t < bucket_of.size())
                bucket_of[(size_t)t] = (int8_t)b;
        }
    }
}

weight_t SpeciesTree::F_bucket_topology(const std::vector<Tree*>& input,
                         const std::array<std::vector<index_t>,4>& buckets,
                         index_t num_taxa,
                         const std::array<int,4>& roles)
{
    std::vector<int8_t> bucket_of;
    build_bucket_of(bucket_of, buckets, num_taxa);

    const std::size_t Sz[4] = {
        buckets[roles[0]].size(),
        buckets[roles[1]].size(),
        buckets[roles[2]].size(),
        buckets[roles[3]].size()
    };

    weight_t r = 0.0;

    std::vector<Node*> order;
    std::vector<std::array<std::size_t,4>> S;

    for (Tree* tree : input) {
        Node* root = tree->get_root();
        postorder_nodes(root, order);

        S.clear();

        for (Node* u : order) {
            if (u->children.size() == 0) {
                // leaf
                index_t t = u->index;
                std::array<std::size_t,4> c = {0,0,0,0};

                int8_t b = ((size_t)t < bucket_of.size()) ? bucket_of[(size_t)t] : (int8_t)-1;
                if (b == roles[0]) c[0] = 1;
                else if (b == roles[1]) c[1] = 1;
                else if (b == roles[2]) c[2] = 1;
                else if (b == roles[3]) c[3] = 1;

                S.push_back(c);
            } else {
                // binary assumed

                auto C2 = S.back(); S.pop_back();
                auto C1 = S.back(); S.pop_back();

                std::array<std::size_t,4> sum = {
                    C1[0] + C2[0],
                    C1[1] + C2[1],
                    C1[2] + C2[2],
                    C1[3] + C2[3]
                };

                std::array<std::size_t,4> rest = {
                    Sz[0] - sum[0],
                    Sz[1] - sum[1],
                    Sz[2] - sum[2],
                    Sz[3] - sum[3]
                };

                std::array<std::array<std::size_t,4>,3> C = {C1, C2, rest};
                r += I_of_C(C);

                S.push_back(sum);
            }
        }
    }

    if (input.empty()) return 0.0;

    return (r * 0.5) / (weight_t)input.size();
}

std::array<weight_t,3> SpeciesTree::freq_three_toplogies(const std::vector<Tree*>& input, Node* blob_node,
                             const std::array<index_t,4>& bucket_ids,
                                       index_t num_taxa)
{   
    index_t temp[4];
    for(index_t i = 0; i < 4; i++) {
        temp[i] = bucket_ids[i];
    }

    quartet_t q = join(temp);
    if (qCFs_average_cache.find(q) != qCFs_average_cache.end()) {
        return qCFs_average_cache[q];
    }

    
    const std::array<int,4> Q1 = {0,1,2,3}; // (0,1)|(2,3)
    const std::array<int,4> Q2 = {0,2,1,3}; // (0,2)|(1,3)
    const std::array<int,4> Q3 = {0,3,1,2}; // (0,3)|(1,2)
    
    std::array<std::vector<index_t>,4> buckets = {blob_node->multi_partitions[bucket_ids[0]],
                                                   blob_node->multi_partitions[bucket_ids[1]],
                                                   blob_node->multi_partitions[bucket_ids[2]],
                                                   blob_node->multi_partitions[bucket_ids[3]]};

    weight_t f1 = F_bucket_topology(input, buckets, num_taxa, Q1);
    weight_t f2 = F_bucket_topology(input, buckets, num_taxa, Q2);
    weight_t f3 = F_bucket_topology(input, buckets, num_taxa, Q3);

    weight_t m = f1 + f2 + f3;
    if (m <= 0.0) return {0.0, 0.0, 0.0};
    qCFs_average_cache[q] = {f1/m, f2/m, f3/m};
    return {f1/m, f2/m, f3/m};
}



// search quard with number of iterations 


weight_t SpeciesTree::search_quard_heuristic(std::vector<Tree *> &input,
                                            std::vector<std::vector<index_t>> &quad,
                                            unsigned long int iter_limit,
                                            index_t *minimizer) {
    index_t indices[4];
    weight_t min = -1;
    size_t count = 0;

    // basic guard
    if (quad.size() != 4) return (weight_t)-1;
    for (int d = 0; d < 4; ++d) if (quad[d].empty()) return (weight_t)-1;

    while (count < iter_limit) {
        // random start: one from each bucket
        indices[0] = quad[0][rand() % quad[0].size()];
        indices[1] = quad[1][rand() % quad[1].size()];
        indices[2] = quad[2][rand() % quad[2].size()];
        indices[3] = quad[3][rand() % quad[3].size()];

        weight_t old_min = min;

        count += neighbor_search_quard(input, quad, indices, &min);

        if (old_min < 0 || min < old_min) {
            minimizer[0] = indices[0];
            minimizer[1] = indices[1];
            minimizer[2] = indices[2];
            minimizer[3] = indices[3];
        }
    }

    return min;
}

weight_t SpeciesTree::search_quard_heuristic(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table,
                                            std::vector<std::vector<index_t>> &quad,
                                            unsigned long int iter_limit,
                                            index_t *minimizer) {
    index_t indices[4];
    weight_t min = -1;
    size_t count = 0;

    // basic guard
    if (quad.size() != 4) return (weight_t)-1;
    for (int d = 0; d < 4; ++d) if (quad[d].empty()) return (weight_t)-1;

    while (count < iter_limit) {
        // random start: one from each bucket
        indices[0] = quad[0][rand() % quad[0].size()];
        indices[1] = quad[1][rand() % quad[1].size()];
        indices[2] = quad[2][rand() % quad[2].size()];
        indices[3] = quad[3][rand() % quad[3].size()];

        weight_t old_min = min;

        count += neighbor_search_quard(qCFs_table, quad, indices, &min);

        if (old_min < 0 || min < old_min) {
            minimizer[0] = indices[0];
            minimizer[1] = indices[1];
            minimizer[2] = indices[2];
            minimizer[3] = indices[3];
        }
    }

    return min;
}



size_t SpeciesTree::neighbor_search_quard(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table,
                                         std::vector<std::vector<index_t>> &quad,
                                         index_t *current,
                                         weight_t *min) {
    auto [current_score, _] = get_pvalue_and_qCFs(qCFs_table,current);
    size_t k = 1;

    while (true) {
        index_t best_next[4] = { current[0], current[1], current[2], current[3] };
        weight_t best_next_score = current_score;

        // Try changing one coordinate at a time (within its bucket)
        for (int d = 0; d < 4; ++d) {
            index_t temp[4] = { current[0], current[1], current[2], current[3] };

            for (size_t i = 0; i < quad[d].size(); ++i) {
                index_t new_index = quad[d][i];
                if (new_index == current[d]) continue;

                temp[d] = new_index;
                auto [temp_score, _] = get_pvalue_and_qCFs(qCFs_table, temp);
                ++k;

                if (temp_score >= 0 && temp_score < best_next_score) {
                    best_next_score = temp_score;
                    best_next[0] = temp[0];
                    best_next[1] = temp[1];
                    best_next[2] = temp[2];
                    best_next[3] = temp[3];
                }

                temp[d] = current[d]; 
            }
        }

        // no improving neighbor
        if (best_next_score >= current_score) break;

        // take best move
        current_score = best_next_score;
        current[0] = best_next[0];
        current[1] = best_next[1];
        current[2] = best_next[2];
        current[3] = best_next[3];
    }

    if (*min < 0 || current_score < *min) *min = current_score;
    return k;
}



size_t SpeciesTree::neighbor_search_quard(std::vector<Tree *> &input,
                                         std::vector<std::vector<index_t>> &quad,
                                         index_t *current,
                                         weight_t *min) {
    auto [current_score, _] = get_pvalue_and_qCFs(input, current);
    size_t k = 1;

    while (true) {
        index_t best_next[4] = { current[0], current[1], current[2], current[3] };
        weight_t best_next_score = current_score;

        // Try changing one coordinate at a time (within its bucket)
        for (int d = 0; d < 4; ++d) {
            index_t temp[4] = { current[0], current[1], current[2], current[3] };

            for (size_t i = 0; i < quad[d].size(); ++i) {
                index_t new_index = quad[d][i];
                if (new_index == current[d]) continue;

                temp[d] = new_index;
                auto [temp_score, _] = get_pvalue_and_qCFs(input, temp);
                ++k;

                if (temp_score >= 0 && temp_score < best_next_score) {
                    best_next_score = temp_score;
                    best_next[0] = temp[0];
                    best_next[1] = temp[1];
                    best_next[2] = temp[2];
                    best_next[3] = temp[3];
                }

                temp[d] = current[d]; 
            }
        }

        // no improving neighbor
        if (best_next_score >= current_score) break;

        // take best move
        current_score = best_next_score;
        current[0] = best_next[0];
        current[1] = best_next[1];
        current[2] = best_next[2];
        current[3] = best_next[3];
    }

    if (*min < 0 || current_score < *min) *min = current_score;
    return k;
}





weight_t SpeciesTree::search_3f1a(std::vector<Tree *> &input, std::tuple<std::vector<Node *>, std::vector<Node *>, std::vector<Node *>, std::vector<Node *>> *quad, index_t* minimizer) {
    index_t i[4] = {0, 0, 0, 0};
    weight_t min = -1;
    
    auto& t = *quad;

    for (index_t alter = 0; alter < 4; alter ++) {
        
        const std::vector<Node*>* current_vec = nullptr;
        
        const std::vector<size_t> partitions_size = {std::get<0>(t).size(), std::get<1>(t).size(), std::get<2>(t).size(), std::get<3>(t).size()};
        

        switch (alter) {
            case 0: current_vec = &std::get<0>(t); break;
            case 1: current_vec = &std::get<1>(t); break;
            case 2: current_vec = &std::get<2>(t); break;
            case 3: current_vec = &std::get<3>(t); break;
        }
        
        index_t temp[4];
        
        for (int j = 0; j < 4; j++) {temp[j] = i[j];}

        for (temp[alter] = 0; temp[alter] < current_vec->size(); temp[alter] ++) {
            index_t cur_quart[4];
            cur_quart[0] = std::get<0>(t)[temp[0]]->index;
            cur_quart[1] = std::get<1>(t)[temp[1]]->index;
            cur_quart[2] = std::get<2>(t)[temp[2]]->index;
            cur_quart[3] = std::get<3>(t)[temp[3]]->index;
            // if it is not the altered partition, we use random index from the partition
            
            // for (int j = 0; j < 4; j++) {
            //     if (j == alter) {
            //         switch (j) {
            //             case 0: cur_quart[0] = std::get<0>(t)[temp[0]]->index; break;
            //             case 1: cur_quart[1] = std::get<1>(t)[temp[1]]->index; break;
            //             case 2: cur_quart[2] = std::get<2>(t)[temp[2]]->index; break;
            //             case 3: cur_quart[3] = std::get<3>(t)[temp[3]]->index; break;
            //         }
            //     } else {
            //         index_t rand_index = rand() % partitions_size[j];
            //         switch (j) {
            //             case 0: cur_quart[0] = std::get<0>(t)[rand_index]->index; break;
            //             case 1: cur_quart[1] = std::get<1>(t)[rand_index]->index; break;
            //             case 2: cur_quart[2] = std::get<2>(t)[rand_index]->index; break;
            //             case 3: cur_quart[3] = std::get<3>(t)[rand_index]->index; break;
            //         }
            //     }
            // }

            weight_t score = get_pvalue(input, cur_quart);
            // auto [score, qcfs] = get_pvalue_and_qCFs(input, cur_quart);
            
            
            if (min < 0 || score < min) {
                i[alter] = temp[alter];
                minimizer[0] = cur_quart[0];
                minimizer[1] = cur_quart[1];
                minimizer[2] = cur_quart[2];
                minimizer[3] = cur_quart[3];
                min = score;
            }
        }
    }
    return min;
}



bool SpeciesTree::is_match_with_split(const std::array<weight_t,3>& qcfs, index_t node_a1_id, index_t node_a2_id, index_t *indices) {
    // determine the split of the quartet based on qcf
    index_t temp[4];
    for (int i = 0; i < 4; ++i) temp[i] = indices[i];
    std::sort(temp, temp + 4);
    std::array<std::array<index_t, 4>, 2> top_2_frequent_topologies = computed_displayed_quartet_toplogy(temp, qcfs);
    std::array<index_t, 4> most_frequent_toplogy = top_2_frequent_topologies[0];
    if (most_frequent_toplogy[0] == node_a1_id && most_frequent_toplogy[1] == node_a2_id) {
        return true;
    } else if (most_frequent_toplogy[0] == node_a2_id && most_frequent_toplogy[1] == node_a1_id) {
        return true;
    } else if (most_frequent_toplogy[2] == node_a1_id && most_frequent_toplogy[3] == node_a2_id) {
        return true;
    } else if (most_frequent_toplogy[2] == node_a2_id && most_frequent_toplogy[3] == node_a1_id) {
        return true;
    } else {
        return false;
    }
}


weight_t SpeciesTree::search_2f2a(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B, index_t* minimizer, size_t &split_match_count, size_t &split_mismatch_count) {
    index_t i[4];
    weight_t min = -1;
    i[0] = 0; i[2] = 0;
        for (i[1] = 1; i[1] < A.size(); i[1] ++) {
                for (i[3] = 1; i[3] < B.size(); i[3] ++) {
                    index_t temp[4];
                    temp[0] = A[i[0]]->index;
                    temp[1] = A[i[1]]->index;
                    temp[2] = B[i[2]]->index;
                    temp[3] = B[i[3]]->index;
                    auto [score, qcfs] = get_pvalue_and_qCFs(input, temp);
                    if (is_match_with_split(qcfs, A[i[0]]->index, A[i[1]]->index, temp)) {
                        split_match_count++;
                    } else {
                        split_mismatch_count++;
                    }
                    
                    if (min < 0 || score < min) {
                        min = score;
                        minimizer[0] = temp[0]; minimizer[1] = temp[1]; minimizer[2] = temp[2]; minimizer[3] = temp[3];
                    }
                }
        }
    return min;
}

weight_t SpeciesTree::search_star(std::vector<Tree *> &input, 
                                  std::vector<Node *> &A,
                                  std::vector<Node *> &B) {
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
                    weight_t score = get_pvalue_star(input, temp);
                    if (min < 0 || score > min) {
                        min = score;
                    }
                }
            }
        }
    }
    //std::cout << "naive iter: " << count << std::endl;
    return min;
}

//// fixed the bug in minimizer updates in the heuristic search
weight_t SpeciesTree::search(std::vector<Tree *> &input, 
                             std::vector<Node *> &A,
                             std::vector<Node *> &B,
                             size_t iter_limit,
                             index_t* minimizer) {
    index_t indices[4];
    weight_t min = -1;
    size_t count = 0;

    while (count < iter_limit) {
        indices[0] = A[rand() % A.size()]->index;
        do {indices[1] = A[rand() % A.size()]->index;} while (indices[0] == indices[1]);
        indices[2] = B[rand() % B.size()]->index;
        do {indices[3] = B[rand() % B.size()]->index;} while (indices[2] == indices[3]);
        
        weight_t old_min = min;

        count += neighbor_search(input, A, B, indices, &min);

        if (min < old_min || old_min < 0) {
            minimizer[0] = indices[0]; minimizer[1] = indices[1]; minimizer[2] = indices[2]; minimizer[3] = indices[3];
        }

        // std::cout << i << ' ' << min << std::endl;
    }
    
    //std::cout << "heuristic iter: " << count << std::endl;
    return min;
}

weight_t SpeciesTree::neighbor_search(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B, index_t *current, weight_t *min) {
    weight_t current_f[3];
    weight_t current_score = get_pvalue(input, current);
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
            temp_score = get_pvalue(input, temp); k ++;
            if (next_score < 0 || temp_score < next_score) {
                next_score = temp_score;
                for (index_t j = 0; j < 4; j ++) 
                    next[j] = temp[j];
                next_f[0] = temp_f[0]; next_f[1] = temp_f[1]; next_f[2] = temp_f[2];
            }
            temp[0] = current[0];
            temp[1] = new_index;
            temp_score = get_pvalue(input, temp); k ++;
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
            temp_score = get_pvalue(input, temp); k ++;
            if (next_score < 0 || temp_score < next_score) {
                next_score = temp_score;
                for (index_t j = 0; j < 4; j ++) 
                    next[j] = temp[j];
                next_f[0] = temp_f[0]; next_f[1] = temp_f[1]; next_f[2] = temp_f[2];
            }
            temp[2] = current[2];
            temp[3] = new_index;
            temp_score = get_pvalue(input, temp); k ++;
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
    }
    return k;
}

weight_t SpeciesTree::search_star(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B, size_t iter_limit) {
    index_t indices[4];
    weight_t min = -1;
    size_t count = 0;
    while (count < iter_limit) {
        indices[0] = A[rand() % A.size()]->index;
        do {indices[1] = A[rand() % A.size()]->index;} while (indices[0] == indices[1]);
        indices[2] = B[rand() % B.size()]->index;
        do {indices[3] = B[rand() % B.size()]->index;} while (indices[2] == indices[3]);
        count += neighbor_search_star(input, A, B, indices, &min);
        // std::cout << i << ' ' << min << std::endl;
    }
    //std::cout << "heuristic iter: " << count << std::endl;
    return min;
}

weight_t SpeciesTree::neighbor_search_star(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B, index_t *current, weight_t *min) {
    weight_t current_f[3];
    weight_t current_score = get_pvalue_star(input, current);
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
            temp_score = get_pvalue_star(input, temp); k ++;
            if (next_score < 0 || temp_score > next_score) {
                next_score = temp_score;
                for (index_t j = 0; j < 4; j ++) 
                    next[j] = temp[j];
                next_f[0] = temp_f[0]; next_f[1] = temp_f[1]; next_f[2] = temp_f[2];
            }
            temp[0] = current[0];
            temp[1] = new_index;
            temp_score = get_pvalue_star(input, temp); k ++;
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
            temp_score = get_pvalue_star(input, temp); k ++;
            if (next_score < 0 || temp_score > next_score) {
                next_score = temp_score;
                for (index_t j = 0; j < 4; j ++) 
                    next[j] = temp[j];
                next_f[0] = temp_f[0]; next_f[1] = temp_f[1]; next_f[2] = temp_f[2];
            }
            temp[2] = current[2];
            temp[3] = new_index;
            temp_score = get_pvalue_star(input, temp); k ++;
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
    }
    return k;
}

// the get_quartet will return the index = the index of sibling of the 0 index in indices - 1 in the indices array, 
// here we alway gonna passing the hybrid to the index of 0 in the indices array, 
// thus we can easily identify the sibilibing of hybrid in the displayed toplogy 
std::pair<weight_t, std::array<weight_t, 3>> SpeciesTree::get_pvalue_and_qCFs(std::vector<Tree *> &input, index_t *indices) {
    
    index_t temp[4];
    for (index_t i = 0; i < 4; i ++) 
        temp[i] = indices[i];
    std::sort(temp, temp + 4);
    quartet_t q = join(temp);

    if (pvalues.find(q) == pvalues.end() || qCFs_cache.find(q) == qCFs_cache.end()) {
        weight_t qCF[3] = {0, 0, 0};
        for (Tree *t : input) {
            index_t topology = t->get_quartet(temp);
            if (topology >= 0) qCF[topology] += 1;
        }

        if ((qCF[0] + qCF[1] + qCF[2]) == 0)
            pvalues[q] = 1.0;
        else
            pvalues[q] = pvalue(qCF);
        qCFs_cache[q] = {qCF[0], qCF[1], qCF[2]}; // main the global qCFs map
    }

    return {pvalues[q], qCFs_cache[q]};
}

std::pair<weight_t, std::array<weight_t, 3>> SpeciesTree::get_pvalue_and_qCFs(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, index_t *indices) {
    
    index_t temp[4];
    for (index_t i = 0; i < 4; i ++) 
        temp[i] = indices[i];
    std::sort(temp, temp + 4);
    quartet_t q = join(temp);


    std::array<weight_t, 3> qCFs;
    qCFs = qCFs_table[q];
    weight_t tmp[3] = {qCFs[0], qCFs[1], qCFs[2]};
    pvalues[q] = pvalue(tmp);
    qCFs_cache[q] = qCFs;
    return {pvalues[q], qCFs};
}


std::array<std::array<index_t, 4>, 2> SpeciesTree::computed_displayed_quartet_toplogy(index_t *indices,
                                                const std::array<weight_t,3>& qcf) {
    index_t temp[4];
    for (int i = 0; i < 4; ++i) temp[i] = indices[i];
    std::sort(temp, temp + 4);

    int best1 = 0;
    for (int t = 1; t < 3; ++t) {
        if (qcf[t] > qcf[best1] || (qcf[t] == qcf[best1] && t < best1)) best1 = t;
    }

    int best2 = (best1 == 0 ? 1 : 0);
    for (int t = 0; t < 3; ++t) {
        if (t == best1) continue;
        if (qcf[t] > qcf[best2] || (qcf[t] == qcf[best2] && t < best2)) best2 = t;
    }

    auto topo_to_quartet = [&](int topo) -> std::array<index_t, 4> {
        // pairs: (ret[0],ret[1]) and (ret[2],ret[3])
        switch (topo) {
            case 0: return { temp[0], temp[1], temp[2], temp[3] }; // 01|23
            case 1: return { temp[0], temp[2], temp[1], temp[3] }; // 02|13
            case 2: return { temp[0], temp[3], temp[1], temp[2] }; // 03|12
            default: return { temp[0], temp[1], temp[2], temp[3] };
        }
    };

    return { topo_to_quartet(best1), topo_to_quartet(best2) };
}



std::array<std::array<index_t, 4>, 2> SpeciesTree::computed_displayed_quartet_toplogy(index_t *indices) {
    index_t temp[4];
    for (index_t i = 0; i < 4; i ++) temp[i] = indices[i];
    std::sort(temp, temp + 4);
    quartet_t q = join(temp);

    index_t best1 = 0, best2 = 1;
    for (index_t t = 0; t < 3; ++t) {
        if (qCFs_cache[t] > qCFs_cache[best1] || (qCFs_cache[t] == qCFs_cache[best1] && t < best1)) best1 = t;
    }
    best2 = (best1 == 0 ? 1 : 0);

    for (index_t t = 0; t < 3; ++t) {
        if (t == best1) continue;
        if (qCFs_cache[t] > qCFs_cache[best2] || (qCFs_cache[t] == qCFs_cache[best2] && t < best2)) best2 = t;
    }

    auto topo_to_quartet = [&](int topo) -> std::array<index_t, 4> {
        // Returns {a,b,c,d} meaning (a,b) and (c,d) are the two sibling pairs.
        switch (topo) {
            case 0: return { temp[0], temp[1], temp[2], temp[3] }; // 01|23
            case 1: return { temp[0], temp[2], temp[1], temp[3] }; // 02|13
            case 2: return { temp[0], temp[3], temp[1], temp[2] }; // 03|12
            default: return { temp[0], temp[1], temp[2], temp[3] };
        }
    };

    return {topo_to_quartet(best1), topo_to_quartet(best2)};
}

static inline index_t sibling_in_topology(const std::array<index_t,4> &topo, index_t taxon) {
    // pair1: topo[0] <-> topo[1], pair2: topo[2] <-> topo[3]
    if (topo[0] == taxon) return topo[1];
    if (topo[1] == taxon) return topo[0];
    if (topo[2] == taxon) return topo[3];
    if (topo[3] == taxon) return topo[2];
    std::cout << "Error: taxon " << taxon << " not found in topology" << std::endl;
    return (index_t)-1;
}

std::array<index_t, 2> SpeciesTree::siblings_in_two_best_topologies(const std::array<std::array<index_t,4>,2> &best2,
                                             index_t taxon) {
    std::array<index_t,2> sibs;
    sibs[0] = sibling_in_topology(best2[0], taxon);
    sibs[1] = sibling_in_topology(best2[1], taxon);
    return sibs;
}



weight_t SpeciesTree::get_pvalue(std::vector<Tree *> &input, index_t *indices) {
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

        if ((qCF[0] + qCF[1] + qCF[2]) == 0)
            pvalues[q] = 1.0;
        else
            pvalues[q] = pvalue(qCF);
        qCFs_cache[q] = {qCF[0], qCF[1], qCF[2]}; // main the global qCFs map
    }
    return pvalues[q];
}

weight_t SpeciesTree::get_pvalue_star(std::vector<Tree *> &input, index_t *indices) {
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
        if ((qCF[0] + qCF[1] + qCF[2]) == 0)
            pvalues_star[q] = -1.0;  // set to -1 so wouldn't impact search for max
        else
            pvalues_star[q] = pvalue_star(qCF);
    }
    return pvalues_star[q];
}

index_t Tree::get_quartet(index_t *indices) {
    Node *leaves[4];
    for (index_t i = 0; i < 4; i ++) {
        // add for absent tadxon in gene tree
        if (!index2node[indices[i]]) {
            return -1;
        }
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
#endif  // ENABLE_TOB