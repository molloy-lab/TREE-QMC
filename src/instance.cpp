#include "instance.hpp"

extern bool DEBUG_MODE;

Instance::Instance(int argc, char **argv) {
    DEBUG_MODE = false;
    input_file = "";
    output_file = "";
    stree_file = "";
    root_str = "";
    normal = "2"; execute = "0"; taxa_mode = "0"; weight = "n"; score_mode = "0";
    support_low = 0; support_high = 0;  // intentionally bad to force user to set
    contract = false; threshold = 0.0;
    refine_seed = 12345; cut_seed = 1; trc = 0; iter_limit = 10;
    
    dict = NULL; output = NULL;
    if (parse(argc, argv)) {
        std::cout << help_info;
    }
    else {
        if (input_file == "") {
            std::cout << "input file not found" << std::endl;
            std::cout << "use -h for more information" << std::endl;
        }
        else {
            dict = new Dict;
            if (input_trees() == 0) {
                std::cout << "Input has " << input.size() << " gene trees and " << dict->size() << " taxa.\n";
                dict->update_singletons();
                resolve_trees();
                prepare_trees();
                // for (Tree *t : input) std::cout << t->to_string() << std::endl;
                if (verbose > "0") {
                    subproblem_csv.open(output_file + "_subproblems.csv");
                    subproblem_csv << "ID,PARENT,DEPTH,SIZE,ARTIFICIAL,SUBSET";
                    if (verbose > "1") {
                        subproblem_csv << ",ENTRY,PRUNED,ZERO";
                    }
                    subproblem_csv << std::endl;
                }
            }
        }
    }
}

Instance::~Instance() {
    for (Tree *t : input) delete t;
    delete output;
    delete dict;
}

long long Instance::solve() {
    if (input.size() == 0) {
        return -1;
    }

    srand(cut_seed);
    std::string mode = normal + execute + taxa_mode + weight;
    auto start = std::chrono::high_resolution_clock::now();

    if (stree_file != "") {
        output = new SpeciesTree(stree_file, dict);
    } else {
        output = new SpeciesTree(input, dict, mode, iter_limit);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    return duration.count();
}

SpeciesTree *Instance::get_solution() {
    return output;
}

void Instance::output_solution() {
    if (score_mode == "1") {
        std::cout << "Computing branch info" << std::endl;
        output->annotate(input);
    }

    if (root_str != "")
        output->root_at_clade(root_taxon_set);
    else 
        output->put_back_root();

    if (output_file == "") {
        if (score_mode == "1")
            std::cout << output->to_string_annotated() << std::endl;
        else
            std::cout << output->to_string_basic() << std::endl;
    }
    else {
        std::ofstream fout(output_file);
        if (score_mode == "1") 
            fout << output->to_string_annotated() << std::endl;
        else
            fout << output->to_string_basic() << std::endl;
        fout.close();
    }
}

bool Instance::parse(int argc, char **argv) {
    std::cout << "wTREE-QMC version 1.0.0\nCOMMAND: ";
    for (int j = 0; j < argc; j ++) 
        std::cout << argv[j] << " ";
    std::cout << std::endl;

    int i = 0;
    bool help = false;
    bool found_support = false;

    while (i < argc) {
        std::string opt(argv[i]);
        if (opt == "-h" || opt == "--help") help = true;

        if (opt == "-i" || opt == "--input") input_file = argv[++ i];
        if (opt == "-o" || opt == "--output") output_file = argv[++ i];
        if (opt == "-q" || opt == "--scoring") stree_file = argv[++ i];
        if (opt == "--root") root_str = argv[++ i];

        if (opt == "-r" || opt == "--support_range") {
            found_support = true;
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2d(param, &support_low) || support_low < 0) {
                std::cout << "ERROR: invalid support value: " << param << "." << std::endl;
                return true;
            }
            param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2d(param, &support_high) || support_high < 0 || support_high <= support_low) {
                std::cout << "ERROR: invalid support value: " << param << "." << std::endl;
                return true;
            }
        }
        if (opt == "-n" || opt == "--normalize") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "0" && param != "1" && param != "2") {
                std::cout << "ERROR: invalid normalize parameter: " << param << "." << std::endl;
                return true;
            }
            normal = param;
        }
        if (opt == "-w" || opt == "--weight") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "n" && param != "s" && param != "h" && param != "l" && param != "f") {
                std::cout << "ERROR: invalid weight parameter: " << param << "." << std::endl;
                return true;
            }
            weight = param;
        }
        if (opt == "-x" || opt == "--execution") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "0" && param != "1" && param != "2" && param != "3" && param != "4") {
                std::cout << "ERROR: invalid execution parameter " << param << "." << std::endl;
                return true;
            }
            execute = param;
        }
        if (opt == "-v" || opt == "--verbose") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "0" && param != "1" && param != "2") {
                std::cout << "ERROR: invalid verbose parameter " << param << "." << std::endl;
                return true;
            }
            verbose = param;
        }
        if (opt == "--polyseed") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2ul(param, &refine_seed)) {
                std::cout << "ERROR: invalid polyseed parameter: " << param << "." << std::endl;
                return true;
            }
        }
        if (opt == "--cutseed") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2ul(param, &cut_seed)) {
                std::cout << "ERROR: invalid cutseed parameter: " << param << "." << std::endl;
                return true;
            }
        }
        if (opt == "--truncate") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2ul(param, &trc)) {
                std::cout << "ERROR: invalid parameter: " << param << "." << std::endl;
                return true;
            }
        }
        if (opt == "--iterlimit") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2ul(param, &iter_limit)) {
                std::cout << "ERROR: invalid parameter: " << param << "." << std::endl;
                return true;
            }
        }
        if (opt == "-c" || opt == "--contract") {
            contract = true;
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2d(param, &threshold)) {
                std::cout << "ERROR: invalid parameter: " << param << "." << std::endl;
                return true;
            }
        }
        if (opt == "--shared") taxa_mode = "1";
        if (opt == "--scored") score_mode = "1";
        i ++;
    }
    std::cout << "input file: " << input_file << std::endl;
    std::cout << "output file: " << (output_file == "" ? "std" : output_file) << std::endl;
    if (stree_file != "")
        std::cout << "score tree file: " << stree_file << std::endl;


    if (root_str != "") {
        // Extract taxa for rooting
        std::cout << "target root placement: " << root_str << std::endl;

        std::vector<int> splits;
        splits.push_back(-1);
        for (int j = 0; j < root_str.size(); j++)
            if (root_str[j] == ',') splits.push_back(j);
        splits.push_back(root_str.size() + 1);

        for (int j = 0; j < splits.size() - 1; j++)
            root_taxon_set.insert(root_str.substr(splits[j]+1, splits[j+1] - splits[j] - 1));

        //for (auto& str : root_taxon_set) std::cout << str << " ";
        //std::cout << std::endl;
    }

    if (contract) {
        std::cout << "contract support threshold: " << (double)threshold << std::endl;
        weight = "n";  // override
    }

    if (weight == "s") {
        std::cout << "weighting mode: support only" << std::endl;
    } else if (weight == "h") {
        std::cout << "weighting mode: hybrid" << std::endl;
    } else if (weight == "l") {
        std::cout << "weighting mode: length only" << std::endl;
    } else if (weight == "f") {
        std::cout << "weighting mode: none fast" << std::endl;
    } else {
        std::cout << "weighting mode: none" << std::endl;
    }

    if (weight == "s" || weight == "h") {
        if (!found_support) {
            std::cout << "ERROR: Must specify branch support range with option -r!" << std::endl;
            return 1;
        }
        std::cout << "support min: " << (double)support_low << ", max: " << (double)support_high << std::endl;
    }

    std::cout << "normalization scheme: n" + normal << std::endl;

    if (execute == "0") {
        std::cout << "execution mode: efficient" << std::endl;
    }
    else if (execute == "1") {
        std::cout << "execution mode: brute force" << std::endl;
    }
    else if (execute == "2") {
        std::cout << "execution mode: compute weighted quartets, then exit" << std::endl;
        std::cout << "quartets saved in: " << output_file << "_quartets.txt" << std::endl;
        quartets_txt.open(output_file + "_quartets.txt");
    }
    else if (execute == "3") {
        std::cout << "execution mode: compute good and bad edges, then exit" << std::endl;
        std::cout << "good edges saved in: " << output_file << "_good_edges.txt" << std::endl;
        std::cout << "bad edges saved in: " << output_file << "_bad_edges.txt" << std::endl;
        good_edges_txt.open(output_file + "_good_edges.txt");
        bad_edges_txt.open(output_file + "_bad_edges.txt");
    }
    else {
        DEBUG_MODE = true;
        execute = "0";
        std::cout << "execution mode: efficient with brute force validation" << std::endl;
    }
    std::cout << "random seed for refinement: " << refine_seed << std::endl;
    std::cout << "random seed for max-cut: " << cut_seed << std::endl;
    std::cout << "max-cut heuristic iteration limit: " << iter_limit << std::endl;

    if (taxa_mode == "1") 
        std::cout << "WARNING: importance values are shared across gene trees!" << std::endl;
    if (trc > 0) 
        std::cout << "WARNING: input tree set is truncated!" << std::endl;

    return help;
}

std::string Instance::get_execution_mode() {
    return execute;
}

std::size_t Instance::input_trees() {
    std::ifstream fin(input_file);
    if (! fin.is_open()) {
        std::cout << "ERROR: input file " << input_file << " does not exist!" << std::endl;
        return 1;
    }
    std::string newick;
    index_t i = 0;
    while (std::getline(fin, newick)) {
        if (newick.find(";") == std::string::npos) break;
        Tree *t = new Tree(newick, dict);
        if (t->size() > 3)
            input.push_back(t);
        else
            std::cout << "Gene tree " << i << " has fewer than 4 species, ignoring." << std::endl;
        if (++ i == trc) break;
    }
    fin.close();
    return 0;
}

void Instance::resolve_trees() {
    srand(refine_seed);
    size_t total = 0;
    for (Tree *t : input)
        total += t->resolve();

    std::cout << "Found " << total << " polytomies." << std::endl;

    if (total > 0 && weight == "f") {
        std::cout << "WARNING: polytomies were refined arbitrarily!" << std::endl;
    }

    if (total != 0) {
        std::ofstream fout(input_file + ".refined");
        for (Tree *t : input)
            fout << t->to_string() << std::endl;
        fout.close();
    }
}

void Instance::prepare_trees() {
    for (Tree *t : input) 
        t->prepare(weight, support_low, support_high, contract, threshold);
}
