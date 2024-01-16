#include "instance.hpp"

extern bool DEBUG_MODE;

Instance::Instance(int argc, char **argv) {
    DEBUG_MODE = false;

    input_file = "";
    output_file = "";
    mapping_file = "";
    stree_file = "";
    table_file = "";
    root_str = "";

    normal_mode = "2";   // use best algorithm for normalizing based on artificial taxa
    execute_mode = "0";  // use fast algorithm
    taxa_mode = "0";     // don't use shared taxon data structure
    weight_mode = "n";   // don't use weighting
    score_mode = "0";    // don't score final tree
    data_mode = "t";     // input data are trees i.e. newick strings
    brln_mode = "g";     // estimate branch lengths under MSC for gene trees
    contract = false;
    char2tree = false;
    rootonly = false;

    support_low = 0;
    support_high = 1.0;
    support_default = 1.0;
    support_threshold = 0.0;

    refine_seed = 12345;
    cut_seed = 1;
    iter_limit = 10;

    dict = NULL;
    output = NULL;

    int outcome = parse(argc, argv);

    if (outcome == 1) {
        std::cout << help_info;
        exit(0);
    }

    if (outcome == 2) {
        std::cout << "use -h for more information" << std::endl;
        exit(1);
    }

    if (outcome > 2) {
        exit(1);
    }

    std::cout << "LOGGING:" << std::endl;

    // Read mapping file
    if (mapping_file != "") prepare_indiv2taxon_map();

    // Process input

    // First try to figure out data type and whether the file exists...

    dict = new Dict;
    if (data_mode == "t") input_trees();
    else input_matrix();

    if (input.size() == 0) {
        std::cout << "\nERROR: Nothing read from input" << std::endl;
        exit(1);
    }

    dict->update_singletons();

    if (char2tree) {
        std::cout << "Writing characters as trees" << std::endl;
        //exit(0);
        std::ofstream fout(output_file);
        if (fout.fail()) {
            for (Tree *t : input) std::cout << t->to_string_basic() << std::endl;
        }
        else {
            for (Tree *t : input) fout << t->to_string_basic() << std::endl;
            fout.close();
        }
        exit(0);
    }

    // TODO: Add suppress uniforcations!
    refine_trees();
    prepare_trees();

    if (verbose > "0") {
        subproblem_csv.open(output_file + "_subproblems.csv");
        subproblem_csv << "ID,PARENT,DEPTH,SIZE,ARTIFICIAL,SUBSET";
        if (verbose > "1") {
            subproblem_csv << ",ENTRY,PRUNED,ZERO";
        }
        subproblem_csv << std::endl;
    }

    // DONE SETTING UP!
}

Instance::~Instance() {
    for (Tree *t : input) delete t;
    delete output;
    delete dict;
}

long long Instance::solve() {
    srand(cut_seed);

    std::string mode = normal_mode + execute_mode + taxa_mode + weight_mode;

    auto start = std::chrono::high_resolution_clock::now();

    if (stree_file != "") {
        output = new SpeciesTree(stree_file, dict);
    } else {
        output = new SpeciesTree(input, dict, mode, iter_limit, output_file);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    return duration.count();
}

SpeciesTree *Instance::get_solution() {
    return output;
}

void Instance::output_solution() {
    std::cout << "Printing species tree" << std::endl;
    std::cout << output->to_string_basic() << std::endl;

    if (root_str != "") {
        std::cout << "Rooting species tree" << std::endl;
        output->root_at_clade(outgroup_taxon_set);
        std::cout << output->to_string_basic() << std::endl;  // want to write for root only
    }
    /*else {
        output->put_back_root();
    }*/

    if (execute_mode == "2" || execute_mode == "3") return;

    if (score_mode == "1") {
        std::cout << "Computing branch info for species tree" << std::endl;

        std::string qfreq_mode = "w";
        if (weight_mode == "n" || weight_mode == "f")
            qfreq_mode = "n";
        output->annotate(input, qfreq_mode);
    }

    if (table_file != "") {
        std::cout << "Writing table" << std::endl;
        std::ofstream fout(table_file);
        if (fout.fail()) {
            std::cout << "  WARNING: Unable to write to " << table_file << std::endl; 
        } else {
            output->write_table(fout, brln_mode);
            fout.close();
        }
    }

    std::cout << "Writing species tree" << std::endl;
    std::ifstream fin(output_file);
    if (!fin.fail()) {
        std::cout << "  WARNING: " << output_file << " already exists, writing to stdout" << std::endl;
        output_file = "";
        fin.close();
    }

    if (output_file != "") {
        std::ofstream fout(output_file);
        if (!fout.fail()) {
            if (rootonly) {
                fout << output->to_string() << std::endl;
	    } else {
                fout << output->to_string_annotated(brln_mode) << std::endl;
	    }
	    fout.close();
            return;
        }
        std::cout << "  WARNING: Unable to write to " << output_file << ", writing to stdout" << std::endl;
    }

    if (rootonly)
        std::cout << output->to_string() << std::endl;
    else
        std::cout << output->to_string_annotated(brln_mode) << std::endl;
}


int Instance::parse(int argc, char **argv) {
    std::cout << "TREE-QMC version 2.0.0" << std::endl;

    std::cout << "COMMAND: ";
    for (int j = 0; j < argc; j ++) 
        std::cout << argv[j] << " ";
    std::cout << std::endl;

    int nweightparam = 0;
    int nminparam = 0;
    int nmaxparam = 0;
    int ndefaultparam = 0;
    bool help = false;

    index_t i = 1;
    while (i < argc) {
        std::string opt(argv[i]);
        if (opt == "-h" || opt == "--help") {
            return 1;
        }
        else if (opt == "-i" || opt == "--input") {
            if (i < argc - 1) input_file = argv[++ i];
        }
        else if (opt == "--chars") {
            data_mode = "c";  // input data are multi-state characters
            brln_mode = "n";  // don't estimate branch lengths!
        }
        else if (opt == "--bp") {
            data_mode = "b";  // input data are 2-state characters
            brln_mode = "b";  // estimate branch lengths under nWF+IS model
        }
        else if (opt == "-a" || opt == "--mapping") {
            if (i < argc - 1) mapping_file = argv[++ i];
        }
        else if (opt == "-u" || opt == "--support") {
            score_mode = "1";
        }
        else if (opt == "-q" || opt == "--supportonly") {
            score_mode = "1";
            if (i < argc - 1) stree_file = argv[++ i];
        }
        else if (opt == "-r" || opt == "--rootonly") {
            rootonly = true;
            if (i < argc - 1) stree_file = argv[++ i];
        }
        else if (opt == "-o" || opt == "--output") {
            if (i < argc - 1) output_file = argv[++ i];
        }
        else if (opt == "--root") {
            if (i < argc - 1) root_str = argv[++ i];
        }
        else if (opt == "--char2tree") {
            char2tree = true;
        }
        else if (opt == "--writetable") {
            if (i < argc - 1) table_file = argv[++ i];
        }

        // Handle weighting options
        else if (opt == "--hybrid") {
            weight_mode = 'h';
            nweightparam += 1;
        }
        else if (opt == "--fast") {
            weight_mode = 'f';
            nweightparam += 1;
        }
        else if (opt == "-w" || opt == "--weight") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param == "1" || param == "h") weight_mode = "h";
            else if (param == "2" || param == "s") weight_mode = "s";
            else if (param == "3" || param == "l") weight_mode = "l";
            else if (param == "4" || param == "n") weight_mode = "n";
            else if (param == "5" || param == "f") weight_mode = "f";
            else {
                std::cout << "\nERROR: invalid (weight) mode parameter: " << param << std::endl;
                return 2;
            }
            nweightparam += 1;
        }

        // Handle support range options
        else if (opt == "-B" || opt == "--bayes") {
            support_default = 0.333;
            support_low = 0.333;
            support_high = 1.0;
            ndefaultparam += 1;
            nminparam += 1;
            nmaxparam += 1;
        }
        else if (opt == "-L" || opt == "--lrt") {
            support_default = 0.0;
            support_low = 0.0;
            support_high = 1.0;
            ndefaultparam += 1;
            nminparam += 1;
            nmaxparam += 1;
        }
        else if (opt == "-S" || opt == "--bootstrap") {
            support_default = 0.0;
            support_low = 0.0;
            support_high = 100.0;
            ndefaultparam += 1;
            nminparam += 1;
            nmaxparam += 1;
        }
        else if (opt == "-d" || opt == "--default") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2d(param, &support_default)) {
                std::cout << "\nERROR: invalid default support value: " << param << std::endl;
                return 2;
            }
            ndefaultparam += 1;
        }
        else if (opt == "-n" || opt == "--min") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2d(param, &support_low)) {
                std::cout << "\nERROR: invalid min support value: " << param << std::endl;
                return 2;
            }
            nminparam += 1;
        }
        else if (opt == "-x" || opt == "--max") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2d(param, &support_high)) {
                std::cout << "\nERROR: invalid max support value: " << param << std::endl;
                return 2;
            }
            nmaxparam += 1;
        }

        // Handle contraction option
        else if (opt == "--contract") {
            contract = true;
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2d(param, &support_threshold)) {
                std::cout << "\nERROR: invalid parameter: " << param << std::endl;
                return 2;
            }
        }

        // Handle artificial taxon normalization options
        else if (opt == "--norm_atax") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "0" && param != "1" && param != "2") {
                std::cout << "\nERROR: invalid graph normalization parameter: " << param << std::endl;
                return 2;
            }
            normal_mode = param;
        }
        else if (opt == "--shared") {
            taxa_mode = "1";
        }

        // Handle other options
        else if (opt == "-v" || opt == "--verbose") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "0" && param != "1" && param != "2") {
                std::cout << "\nERROR: invalid verbose parameter " << param << std::endl;
                return 2;
            }
            verbose = param;
        }
        else if (opt == "-e" || opt == "--execute") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "0" && param != "1" && param != "2" && param != "3" && param != "4") {
                std::cout << "\nERROR: invalid execution parameter " << param << std::endl;
                return 2;
            }
            execute_mode = param;
        }
        else if (opt == "--polyseed") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2ul(param, &refine_seed)) {
                std::cout << "\nERROR: invalid polyseed parameter: " << param << std::endl;
                return 2;
            }
        }
        else if (opt == "--cutseed") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2ul(param, &cut_seed)) {
                std::cout << "\nERROR: invalid cutseed parameter: " << param << std::endl;
                return 2;
            }
        }
        else if (opt == "--iterlimit") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2ul(param, &iter_limit)) {
                std::cout << "\nERROR: invalid parameter: " << param << std::endl;
                return 2;
            }
        }
        else {
            std::cout << "WARNING: Ignoring unrecognized option: " << opt << std::endl;
        }

        i ++;
    }

    // Report settings and check that they make sense
    std::cout << std::endl << "SETTINGS:" << std::endl;

    // Input files
    if (input_file == "") {
        std::cout << "input file not found" << std::endl;
        return 2;   
    }
    std::cout << "input file: " << input_file << std::endl;
    if (mapping_file != "") std::cout << "mapping file: " << mapping_file << std::endl;
    if (stree_file != "") std::cout << "species tree file: " << stree_file << std::endl;
    if (table_file != "") std::cout << "table file: " << table_file << std::endl;

    // Output file
    if (output_file == "")
         std::cout << "output file: std" << std::endl;
    else
        std::cout << "output file: " << output_file << std::endl;

    // Get outgroup taxon set for rooting
    if (root_str != "") {
        prepare_root_taxa();
        std::cout << "target root placement:";
        for (auto &taxon : outgroup_taxon_set) std::cout << " " << taxon;
        std::cout << std::endl;
    }

    if (data_mode == "t") {
        // Process weighting options
        if (nweightparam > 1) {
            std::cout << "\nERROR: Conflicting weighting options specified" << std::endl;
            return 2;
        }
        if (weight_mode == "s")
            std::cout << "weighting mode: support only" << std::endl;
        else if (weight_mode == "h")
            std::cout << "weighting mode: hybrid" << std::endl;
        else if (weight_mode == "l")
            std::cout << "weighting mode: length only" << std::endl;
        else if (weight_mode == "f")
            std::cout << "weighting mode: none (fast)" << std::endl;
        else if (weight_mode == "n")
            std::cout << "weighting mode: none" << std::endl;
        else {
            std::cout << "\nERROR: Unrecognized weight option!" << std::endl;
            return 2;
        }
        if (weight_mode != "h")
            std::cout << "  WARNING: --hybrid option is recommended" << std::endl;

        // Process support branch options
        if (weight_mode == "s" || weight_mode == "h" || contract) {
            // Check support options make sense
            if (nminparam  == 0 || nmaxparam == 0 || ndefaultparam == 0) {
                std::cout << "\nERROR: Must specify min, max, and default support values or use preset option" << std::endl;
                return 2;
            }
            else if (nminparam > 1 || nmaxparam > 1 || ndefaultparam > 1) {
                std::cout << "\nERROR: Multiple support parameters specified" << std::endl;
                return 2;
            }

            std::cout << "support min: " << (double)support_low << std::endl;
            std::cout << "support max: " << (double)support_high << std::endl;
            std::cout << "support default: " << (double)support_default << std::endl;

            // Check support parameters make sense 
            if (support_low >= support_high || support_high <= support_low) {
                std::cout << "\nERROR: Conflicting support values specified" << std::endl;
                return 3;
            }
            if (support_low < 0) {
                std::cout << "\nERROR: Support min value must be non-negative" << std::endl;
                return 3;
            }
            if (support_default < support_low || support_default > support_high) {
                std::cout << "\nERROR: Default support value is outside of range" << std::endl;
                return 3;
            }
        } else {
            // Ignore any support values that were specified and set back to defaults
            if (nminparam > 0 || nmaxparam > 0 || ndefaultparam > 0)
                std::cout << "  WARNING: Running in unweighted mode (no contraction) so ignoring support options" << std::endl;
            support_low = 0.0;
            support_high = 1.0;
            support_default = 1.0;
        }

        // Process contraction options
        if (contract) {
            if (weight_mode == "f") {
                std::cout << "\nERROR: Cannot contract low support branches in --fast mode" << std::endl;
                return 2;
            }

            std::cout << "contract support threshold: " << (double)support_threshold << std::endl;

            if (support_threshold < 0 || support_threshold > 1) {
                std::cout << "\nERROR: Support threshold must be between 0 and 1 (it's applied after support values are mapped to this interval)" << std::endl;
                return 1;
            }
        }
    }
    else {
        if (nweightparam > 1 || nminparam > 0 || nmaxparam > 0 || ndefaultparam > 0) {
            if (data_mode == "b")
                std::cout << "  WARNING: Running in bipartition mode so ignoring any weight or support options" << std::endl;
            else
                std::cout << "  WARNING: Running in character mode so ignoring any weight or support options" << std::endl;
        }
        weight_mode = "n";
        support_low = 0.0;
        support_high = 1.0;
        support_default = 1.0;
    }

    // Process normalization mode
    std::cout << "normalization for artifical taxa mode: n" + normal_mode;
    if (taxa_mode == "1") std::cout << " (shared)";
    std::cout << std::endl;
    if (normal_mode != "2")
        std::cout << "  WARNING: --norm_atax 2 is recommended" << std::endl;

    // Process execution mode
    if (execute_mode == "0") {
        std::cout << "execution mode: efficient" << std::endl;
    }
    else if (execute_mode == "1") {
        std::cout << "execution mode: brute force" << std::endl;
    }
    else if (execute_mode == "2") {
        std::cout << "execution mode: compute weighted quartets, then exit" << std::endl;
        std::cout << "good edges will be saved in: " << output_file << "_quartets.txt" << std::endl;
    }
    else if (execute_mode == "3") {
        std::cout << "execution mode: compute good and bad edges, then exit" << std::endl;
        std::cout << "good edges will be saved in: " << output_file << "_good_edges.txt" << std::endl;
        std::cout << "bad edges will be saved in: " << output_file << "_bad_edges.txt" << std::endl;
    }
    else {
        DEBUG_MODE = true;
        execute_mode = "0";
        std::cout << "execution mode: efficient with brute force validation" << std::endl;
    }

    // Print remaining options
    std::cout << "random seed for refinement: " << refine_seed << std::endl;
    std::cout << "random seed for max-cut: " << cut_seed << std::endl;
    std::cout << "max-cut heuristic iteration limit: " << iter_limit << std::endl;

    std::cout << std::endl;

    return 0;
}

std::string Instance::get_execution_mode() {
    return execute_mode;
}

void Instance::input_trees() {
    std::ifstream fin(input_file);
    if (fin.fail()) {
        std::cout << "\nERROR: Unable to open " << input_file << std::endl;
        exit(1);
    }

    std::cout << "Reading input" << std::endl;
    std::string newick;
    int mintax = INDEX_WIDTH;
    int maxtax = 0;
    index_t i = 1;
    while (std::getline(fin, newick)) {
        // TODO: change to function that checks if newick string is valid, before proceeding
        if (newick.find(";") != std::string::npos) {
            Tree *t = new Tree(newick, dict, indiv2taxon, support_default);
            if (t->size() > maxtax) maxtax = t->size();
            if (t->size() < mintax) mintax = t->size();
            if (t->size() > 3)
                input.push_back(t);
            else
                std::cout << "  WARNING: Input tree on line " << i << " has fewer than 4 species so ignoring" << std::endl;
        }
        i++;
    }

    std::cout << "Found" << std::endl;
    std::cout << "    " << input.size() << " trees\n";
    std::cout << "    " << dict->size() << " taxa\n";

    if (mintax != maxtax && taxa_mode == "1") {
        std::cout << "    some input trees are missing taxa" << std::endl;
        std::cout << "WARNING: --shared option should NOT be used with missing taxa!" << std::endl;
    }

    fin.close();
}

void Instance::input_matrix() {
    std::ifstream fin(input_file);
    if (fin.fail()) {
        std::cout << "\nERROR: Unable to open " << input_file << std::endl;
        exit(1);
    }
    std::string line;
    std::getline(fin, line);
    fin.close();

    if (line[0] == '(') {
        input_trees();
        return;
    }

    CharMat *cmat = new CharMat(input_file);

    std::string newick;
    int mintax = INDEX_WIDTH;
    int maxtax = 0;
    index_t i = 1;
    while (cmat->size() > 0) {
        cmat->pop_newick(newick);
        if (newick != "") {
            //std::cout << newick << std::endl;
            Tree *t = new Tree(newick, dict, indiv2taxon, support_default);
            //std::cout << t->to_string_basic() << std::endl;
            //exit(1);
            if (t->size() > maxtax) maxtax = t->size();
            if (t->size() < mintax) mintax = t->size();
            input.push_back(t);
        }
        i++;
    }


    std::cout << "Found" << std::endl;
    std::cout << "    " << input.size() << " informative characters\n";
    std::cout << "    " << dict->size() << " taxa\n";

    if (mintax != maxtax && taxa_mode == "1") {
        std::cout << "    some input characters are missing taxa" << std::endl;
        std::cout << "WARNING: --shared option should NOT be used with missing taxa!" << std::endl;
    }
}

void Instance::refine_trees() {
    srand(refine_seed);

    std::size_t total = 0;
    for (Tree *t : input)
        total += t->refine();

    if (data_mode == "t")
        std::cout << "    " << total << " polytomies across all trees\n";

    if (total > 0 && weight_mode == "f") {
        std::cout << "  WARNING: polytomies were refined arbitrarily!" << std::endl;

        // Save refined trees for inspection!
        std::ofstream fout(input_file + ".refined");
        if (!fout.fail()) {
            for (Tree *t : input)
                fout << t->to_string() << std::endl;
            fout.close();
        }
    }
}

void Instance::prepare_trees() {
    for (Tree *t : input) 
        t->prepare(weight_mode, support_low, support_high, support_threshold);
}

void Instance::prepare_root_taxa() {
    std::string taxon;

    // Try to read outgroup taxon set from file (1 taxon per line)
    std::ifstream fin(root_str);
    if (!fin.fail()) {
        while (std::getline(fin, taxon)) {
            if (taxon != "")
                outgroup_taxon_set.insert(taxon);
        }
        fin.close();
        return;
    }

    // Read outgroup from comma separated string
    // TODO: use string stream with getline and delim=','
    std::vector<int> splits;

    splits.push_back(-1);
    for (int j = 0; j < root_str.size(); j++)
        if (root_str[j] == ',')
            splits.push_back(j);
    splits.push_back(root_str.size() + 1);

    for (int j = 0; j < splits.size() - 1; j++) {
        taxon = root_str.substr(splits[j]+1, splits[j+1] - splits[j] - 1);
        outgroup_taxon_set.insert(taxon);
    }
}


void Instance::prepare_indiv2taxon_map() {
    std::string taxon, indiv, line, delim;
    std::size_t sep;

    std::ifstream fin(mapping_file);
    if (fin.fail()) {
        std::cout << "\nERROR: Unable to open mapping file " << mapping_file << std::endl;
        exit(1);
    }

    std::cout << "Reading mapping file" << std::endl;

    // Process first line and figure out the delimitator
    std::getline(fin, line);
    delim = " ";
    sep = line.find(delim);
    if (sep == std::string::npos) {
        delim = "\t";
        sep = line.find(delim);
        if (sep == std::string::npos) {
            std::cout << "\nERROR: Unable to process mapping file" << std::endl;
            exit(1);
        }
    }
    if (sep == 0 || sep + 1 == std::string::npos) {
        std::cout << "  WARNING: Unable to process mapping file line " << line << std::endl;
    }
    else {
        indiv = line.substr(0, sep);
        taxon = line.substr(sep + 1, std::string::npos);
        indiv2taxon.insert({indiv, taxon});
    }

    // Process remainder of file
    while (std::getline(fin, line)) {
        sep = line.find(delim);
        if (sep == std::string::npos || sep == 0 || sep + 1 == std::string::npos) {
            std::cout << "  WARNING: Unable to process mapping file line " << line << std::endl;
        }
        else {
            indiv = line.substr(0, sep);
            taxon = line.substr(sep + 1, std::string::npos);
            indiv2taxon.insert({indiv, taxon});
        }
    }

    //for (auto itr = indiv2taxon.begin(); itr != indiv2taxon.end(); ++itr)
    //    std::cout << itr->first << " : " << itr->second << std::endl;

    fin.close();
}
