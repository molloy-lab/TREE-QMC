#include "graph.hpp"
#include <algorithm>
#include <cstdlib>
#include <vector>

#ifdef TREEQMC_USE_OPENMP
#include <omp.h>
#endif

namespace {

int treeqmc_env_int(const char *name, int fallback) {
    const char *env_p = std::getenv(name);
    if (env_p == NULL || *env_p == '\0') return fallback;
    const int value = std::atoi(env_p);
    return value > 0 ? value : fallback;
}

unsigned long long treeqmc_env_ull(const char *name, unsigned long long fallback) {
    const char *env_p = std::getenv(name);
    if (env_p == NULL || *env_p == '\0') return fallback;
    const unsigned long long value = std::strtoull(env_p, nullptr, 10);
    return value > 0 ? value : fallback;
}

#ifdef TREEQMC_USE_OPENMP
int treeqmc_graph_threads(std::size_t tree_count) {
    const int requested = treeqmc_env_int("TREEQMC_NUM_THREADS", 0);
    if (requested > 0) return std::max(1, std::min<int>(requested, static_cast<int>(tree_count)));
    const int max_threads = omp_get_max_threads();
    return std::max(1, std::min<int>(max_threads, std::min<int>(64, static_cast<int>(tree_count))));
}
#endif

}  // namespace

Graph::Graph(std::vector<Tree *> trees, Taxa &subset, std::string weighting) {
    size = subset.size();
    for (index_t i = 0; i < size; i ++) {
        index2index[subset.root_at(i)] = i;
        indices.push_back(subset.root_at(i));
    }
    graph = new weight_t**[2];
    graph[0] = Matrix::new_mat(size);
    graph[1] = Matrix::new_mat(size);
    if (verbose > "1") count[1] = count[2] = count[3] = 0;

    const std::size_t dense_n = static_cast<std::size_t>(size) * static_cast<std::size_t>(size);
    weight_t ***subgraph = new weight_t**[2];
    subgraph[0] = Matrix::new_mat(size);
    subgraph[1] = Matrix::new_mat(size);

    weight_t *graph_good_ptr = graph[0][0];
    weight_t *graph_bad_ptr = graph[1][0];

#ifdef TREEQMC_USE_OPENMP
    const unsigned long long estimated_work =
        static_cast<unsigned long long>(dense_n) *
        static_cast<unsigned long long>(trees.size());
    const unsigned long long parallel_min_work =
        treeqmc_env_ull("TREEQMC_PARALLEL_MIN_WORK", 1ULL);
    const int parallel_threads = treeqmc_graph_threads(trees.size());
    const bool use_parallel =
        parallel_threads > 1 &&
        trees.size() > 1 &&
        verbose <= "1" &&
        estimated_work >= parallel_min_work;

    if (use_parallel) {
#pragma omp parallel num_threads(parallel_threads)
        {
            Taxa local_subset(subset);
            weight_t ***local_subgraph = new weight_t**[2];
            local_subgraph[0] = Matrix::new_mat(size);
            local_subgraph[1] = Matrix::new_mat(size);
            std::vector<weight_t> local_good(dense_n, 0.0);
            std::vector<weight_t> local_bad(dense_n, 0.0);

#pragma omp for schedule(static)
            for (std::size_t tree_idx = 0; tree_idx < trees.size(); ++tree_idx) {
                Tree *tree = trees[tree_idx];
                std::unordered_map<index_t, index_t> &valid = tree->get_indices();
                local_subset.weight_update(valid);
                if (weighting == "f")
                    tree->build_graph_into(local_subset, local_subgraph);
                else
                    tree->build_wgraph_into(local_subset, local_subgraph);

                weight_t *sub_good = local_subgraph[0][0];
                weight_t *sub_bad = local_subgraph[1][0];
                for (std::size_t idx = 0; idx < dense_n; ++idx) {
                    const weight_t good = sub_good[idx];
                    const weight_t bad = sub_bad[idx];
                    if (good > 0 || bad > 0) {
                        local_good[idx] += good;
                        local_bad[idx] += bad;
                    }
                }
            }

#pragma omp critical
            {
                for (std::size_t idx = 0; idx < dense_n; ++idx) {
                    graph_good_ptr[idx] += local_good[idx];
                    graph_bad_ptr[idx] += local_bad[idx];
                }
            }

            Matrix::delete_mat(local_subgraph[0], size);
            Matrix::delete_mat(local_subgraph[1], size);
            delete [] local_subgraph;
        }
    }
    else
#endif
    {
        auto accumulate_subgraph_to_graph = [&]() {
            weight_t *sub_good = subgraph[0][0];
            weight_t *sub_bad = subgraph[1][0];
            for (std::size_t idx = 0; idx < dense_n; ++idx) {
                const weight_t good = sub_good[idx];
                const weight_t bad = sub_bad[idx];
                if (good > 0 || bad > 0) {
                    graph_good_ptr[idx] += good;
                    graph_bad_ptr[idx] += bad;
                }
            }
        };

        for (Tree *tree : trees) {
            std::unordered_map<index_t, index_t> &valid = tree->get_indices();
            subset.weight_update(valid);
            if (weighting == "f")
                tree->build_graph_into(subset, subgraph);
            else
                tree->build_wgraph_into(subset, subgraph);
            accumulate_subgraph_to_graph();
        }
    }

    Matrix::delete_mat(subgraph[0], size);
    Matrix::delete_mat(subgraph[1], size);
    delete [] subgraph;
    /*
    weight_t **temp_graph = Matrix::new_mat(size);
    for (index_t i = 0; i < size; i ++) {
        for (index_t j = 0; j < size; j ++) {
            temp_graph[i][j] = graph[0][i][j] + graph[1][i][j];
        }
    }
    std::cout << Matrix::display_mat(temp_graph, size) << std::endl;
    Matrix::delete_mat(temp_graph, size);
    */
    if (verbose > "1") subproblem_csv << ',' << count[1] << ',' << count[2] << ',' << count[3];
}

void Graph::write_good_edges(Dict *dict) {
    index_t i_, j_;
    weight_t normval, nelem;
    
    // Find average value in good edges   
    /*nelem = size * (size - 1) / 2; // 0 + 1 + 2 + ... + n-1 = n(n-1)/2
    normval = 0.0;
    for (size_t i = 1; i < size-1; i++) {
        for (size_t j = 0; j < i; j++) {
            normval += (graph[0][i][j] / nelem);
        }
    }*/

    // Find minimum value in good edges
    normval = std::numeric_limits<double>::max();
    for (size_t i = 1; i < size-1; i++) {
        for (size_t j = 0; j < i; j++) {
            if (graph[0][i][j] > 1 && graph[0][i][j] < normval)
                normval = graph[0][i][j];
        }
    }

    std:: cout << "Normalizing good edges by " << normval << std::endl;

    good_edges_txt << size << std::endl;
    for (size_t i = 0; i < size; i++) {
        i_ = index2index[i];
        good_edges_txt << dict->index2label(i);
        for (size_t j = 0; j < size; j++) {
            j_ = index2index[j];
            good_edges_txt << " " << graph[0][i_][j_] / normval;
        }
        good_edges_txt << std::endl;
    }
}

void Graph::write_bad_edges(Dict *dict) {
    index_t i_, j_;
    weight_t normval, nelem;
    
    // Find average value in good edges   
    /*nelem = size * (size - 1) / 2; // 0 + 1 + 2 + ... + n-1 = n(n-1)/2
    normval = 0.0;
    for (size_t i = 1; i < size-1; i++) {
        for (size_t j = 0; j < i; j++) {
            normval += (graph[1][i][j] / nelem);
        }
    }*/

    // Find minimum value in good edges
    normval = std::numeric_limits<double>::max();
    for (size_t i = 1; i < size-1; i++) {
        for (size_t j = 0; j < i; j++) {
            if (graph[1][i][j] > 1 && graph[1][i][j] < normval)
                normval = graph[1][i][j];
        }
    }

    std:: cout << "Normalizing bad edges by " << normval << std::endl;

    bad_edges_txt << size << std::endl;
    for (size_t i = 0; i < size; i++) {
        i_ = index2index[i];
        bad_edges_txt << dict->index2label(i);
        for (size_t j = 0; j < size; j++) {
            j_ = index2index[j];
            bad_edges_txt << " " << graph[1][i_][j_] / normval;
        }
        bad_edges_txt << std::endl;
    }
}

Graph::Graph(std::unordered_map<quartet_t, weight_t> &quartets, Taxa &subset) {
    size = subset.size();
    for (index_t i = 0; i < size; i ++) {
        index2index[subset.root_at(i)] = i;
        indices.push_back(subset.root_at(i));
    }
    graph = new weight_t**[2];
    graph[0] = Matrix::new_mat(size);
    graph[1] = Matrix::new_mat(size);
    for (auto elem : quartets) {
        index_t *indices = split(elem.first);
        index_t a = index2index[indices[0]], b = index2index[indices[1]], c = index2index[indices[2]], d = index2index[indices[3]];
        weight_t w = elem.second;
        // bad edges
        graph[1][a][b] += w; graph[1][c][d] += w; graph[1][b][a] += w; graph[1][d][c] += w;

        // good edges
        graph[0][a][c] += w; graph[0][a][d] += w; graph[0][b][c] += w; graph[0][b][d] += w;
        graph[0][c][a] += w; graph[0][d][a] += w; graph[0][c][b] += w; graph[0][d][b] += w;
        delete [] indices;
    }
    /*
    weight_t **temp_graph = Matrix::new_mat(size);
    for (index_t i = 0; i < size; i ++) {
        for (index_t j = 0; j < size; j ++) {
            temp_graph[i][j] = graph[0][i][j] + graph[1][i][j];
        }
    }
    std::cout << Matrix::display_mat(temp_graph, size) << std::endl;
    Matrix::delete_mat(temp_graph, size);
    */
}

Graph::~Graph() {
    Matrix::delete_mat(graph[0], size);
    Matrix::delete_mat(graph[1], size);
    delete[] graph;
}

std::string Graph::to_string() {
    return Matrix::display_mat(graph[0], size) + "\n" + Matrix::display_mat(graph[1], size);
}

weight_t Graph::get_cut(std::vector<index_t> *A, std::vector<index_t> *B, unsigned long int iter_limit) {
    weight_t positive_weight = -1.0;
    std::vector<index_t> a, b;
    weight_t lower = 0.0, upper = 6.0;
    /*
    for (index_t i = 0; i < size; i ++) {
        for (index_t j = i + 1; j < size; j ++) {
            if (graph[1][i][j] == 0) continue;
            weight_t ratio = graph[0][i][j] / graph[1][i][j];
            if (ratio > upper) upper = ratio;
        }
    }
    */
    while (lower + 0.1 < upper) {
        weight_t alpha = (lower + upper) / 2.0;
        a.clear(); b.clear();
        weight_t weight = sdp_cut(alpha, &a, &b, iter_limit);
        if (weight < 0.001 || a.size() <= 1 || b.size() <= 1) {
            upper = alpha;
        }
        else {
            lower = alpha;
            positive_weight = alpha;
            *A = a;
            *B = b;
        }
    }
    /*
    if (A->size() <= 1 || B->size() <= 1) {
        std::cout << Matrix::display_mat(graph[0], size) << std::endl;
        std::cout << Matrix::display_mat(graph[1], size) << std::endl;
    }
    assert(A->size() > 1 && B->size() > 1);
    */
    return positive_weight;
}

weight_t Graph::sdp_cut(weight_t alpha, std::vector<index_t> *A, std::vector<index_t> *B, unsigned long int iter_limit) {
    std::vector<Instance::InstanceTuple> input;
    weight_t avg = 0, num = size * (size - 1) / 2;
    for (index_t i = 0; i < size; i ++) {
        for (index_t j = i + 1; j < size; j ++) {
            weight_t temp = (graph[0][i][j] - alpha * graph[1][i][j]) / num;
            if (temp < 0) temp = - temp;
            avg += temp;
        }
    }
    for (index_t i = 0; i < size; i ++) {
        for (index_t j = i + 1; j < size; j ++) {
            weight_t weight = (graph[0][i][j] - alpha * graph[1][i][j]) / avg;
            input.push_back(Instance::InstanceTuple(std::make_pair(i + 1, j + 1), weight));
        }
    }
    QuartetGraphMaxCutCallback mc(iter_limit);
    MaxCutInstance instance(input, size);
    Burer2002 heuristic(instance, -1, false, &mc);
    MaxCutSimpleSolution solution = heuristic.get_best_solution();
    std::vector<int> cut = solution.get_assignments();
    for (index_t i = 0; i < cut.size(); i ++) {
        if (cut[i] < 0) 
            A->push_back(indices[i]);
        else 
            B->push_back(indices[i]);
    }
    return solution.get_weight();
}

QuartetGraphMaxCutCallback::QuartetGraphMaxCutCallback(unsigned long int iter_limit) {
    this->iter_limit = iter_limit;
}

QuartetGraphMaxCutCallback::~QuartetGraphMaxCutCallback() {

}

bool QuartetGraphMaxCutCallback::Report(const MaxCutSimpleSolution& solution, bool newBest, double runtime) {
    return true;
}

bool QuartetGraphMaxCutCallback::Report(const MaxCutSimpleSolution& solution, bool newBest, double runtime, int iter) {
    assert(iter >= 0);
    return iter < this->iter_limit;
}
