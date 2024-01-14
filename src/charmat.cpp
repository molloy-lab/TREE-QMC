#include "charmat.hpp"

CharMat::CharMat(std::string input_file) {
	std::ifstream fin(input_file);
    if (fin.fail()) {
        std::cout << "\nERROR: Unable to open " << input_file << std::endl;
        exit(1);
    }
    std::string line;
    std::getline(fin, line);
    fin.seekg(0, std::ios::beg);  // roll back to beginning
    

    if (line[0] == '>')
        read_fasta(fin);
    else if (line[0] == '#')
        read_nexus(fin);
    else
        read_phylip(fin);

    fin.close();
}

CharMat::~CharMat() {}

void CharMat::pop_newick(std::string &newick) {
    std::string split;
    std::size_t ntax, ntax_done, ntax_found, ntax_notmiss, test;
    char c;

    split = splits_.back();
    ntax = split.size();

    if (ntax != labels_.size()) {
        std::cout << "ERROR: Unable to process character matrix!" << std::endl;
        exit(1);
    }

    // dictionary maps states to pair
    // first - comma separated list of taxa in that state
    // second - number of taxa in state
    std::unordered_map<char, std::pair<std::string, int>> mymap;
    for (size_t j = 0; j < ntax; j++) {
        c = split[j];

        auto itr = mymap.find(c);
        if (itr == mymap.end()) {
            mymap[c] = std::make_pair(labels_[j], 1);
        }
        else {
            itr->second.first.append("," + labels_[j]);
            itr->second.second++;
            //std::cout << itr->first << ": " << itr->second.first << " (" << itr->second.second << " taxa)" << std::endl;
        }
    }

    // Convert to newick string
    ntax_done = 0;
    ntax_notmiss = 0;
    test = 0;  // checking at least 2 states with 2 or more taxa
    newick = "(";

    // Remove missing taxa
    auto itr = mymap.find('?');
    if (itr != mymap.end()) {
        ntax_done += itr->second.second;
        mymap.erase('?');
    }
    itr = mymap.find('-');
    if (itr != mymap.end()) {
        ntax_done += itr->second.second;
        mymap.erase('-');
    }
    itr = mymap.find('N');
    if (itr != mymap.end()) {
        ntax_done += itr->second.second;
        mymap.erase('N');
    }

    // Process ancestral state
    itr = mymap.find('0');
    if (itr != mymap.end()) {
        newick += itr->second.first;
        ntax_found = itr->second.second;

        if (ntax_found > 1) test++;
        ntax_notmiss += ntax_found;

        ntax_done += ntax_found;
        if (ntax_done != ntax) newick.append(",");

        mymap.erase('0');
    }

    for (itr = mymap.begin(); itr != mymap.end(); ++itr) {
        ntax_found = itr->second.second;
        ntax_notmiss += ntax_found;

        if (ntax_found == 1) {
            newick.append(itr->second.first);
        } else {
            newick.append("(");
            newick.append(itr->second.first);
            newick.append(")");
            test++;
        }

        ntax_done += ntax_found;
        if (ntax_done != ntax) newick.append(",");
    }

    newick.append(");");

    if (test < 2) newick = "";
    if (ntax_notmiss < 4) newick = "";

    //std::cout << split << " vs " << newick << std::endl << std::endl;

    // Pop split
    splits_.pop_back();
}

std::size_t CharMat::size() {
    return splits_.size();
}

void CharMat::read_fasta(std::istream &is) {
    std::cout << "Reading input as fasta file" << std::endl;

    std::string line;
    std::size_t nseq, nchr;
    char c;
    bool nextlab;

    // Confirm fasta format
    std::getline(is, line);
    if (line[0] != '>') {
        std::cout << "ERROR: Unable to read as fasta file" << std::endl;
        exit(1);
    }
    labels_.push_back(line.substr(1));
    nseq = 1;

    // Count characters
    std::getline(is, line);
    if (line[0] == '>') {
        std::cout << "ERROR: Unable to read as fasta file" << std::endl;
        exit(1);
    }
    nchr = line.size();

    // Set up split matrix
    for (int i = 0; i < nchr; i++) {
        splits_.push_back("");
        splits_[i].push_back(line[i]);
    }

    // Read data
    nextlab = true;
    while (std::getline(is, line)) {
        if (line[0] == '>') {
            // Read label
            if (!nextlab) {
                std::cout << "ERROR: Unable to read as fasta file" << std::endl;
                exit(1);
            }
            labels_.push_back(line.substr(1));
            nextlab = false;
            nseq++;
        }
        else {
            // Read sequence data
            if (nextlab) {
                std::cout << "ERROR: Unable to read as fasta file" << std::endl;
                exit(1);
            }
            if (line.length() != nchr) {
                std::cout << "ERROR: Wrong sequence length for " << labels_[nseq-1] << std::endl;
                exit(1);
            }
            for (int i = 0; i < nchr; i++) {
                splits_[i].push_back(line[i]);
            }
            nextlab = true;
        }
    }

    //for (int i = 0; i < splits_.size(); i++)
    //    std::cout << splits_[i] << std::endl;

    //for (int i = 0; i < labels_.size(); i++)
    //    std::cout << labels_[i] << std::endl;
}

void CharMat::read_nexus(std::istream &is) {
    std::cout << "Reading input as nexus file" << std::endl;

    std::string line, word;
    std::size_t sep, pos, spos, epos;
    unsigned long ntax, nchr;
    bool found;

    // Roll forward until you find Dimensions ntax= nchar=
    found = false;
    while (std::getline(is, line)) {
        for (auto & c: line) c = std::tolower(c);
        sep = line.find("dimensions");
        if (sep != std::string::npos) {
            found = true;
            break;
        }
    }
    if (!found) {
        std::cout << "ERROR: Unable to find dimensions block of nexus file"<< std::endl;
        exit(1);
    }

    // Parse number of taxa
    spos = line.find("ntax=");
    if (spos == std::string::npos) {
        std::cout << "ERROR: Unable to find ntax variable in nexus file"<< std::endl;
        exit(1);
    }
    word = line.substr(spos + 5, std::string::npos);
    epos = word.find(" ");
    if (epos == std::string::npos) epos = word.find(";");
    word = word.substr(0, epos);
    if (! s2ul(word, &ntax)) {
        std::cout << "ERROR: Number of taxa is not an integer"<< std::endl;
        exit(1);
    }    

    // Parse number of characters
    spos = line.find("nchar=");
    if (spos == std::string::npos) {
        std::cout << "ERROR: Unable to find nchar variable in nexus file"<< std::endl;
        exit(1);
    }
    word = line.substr(spos + 6, std::string::npos);
    epos = word.find(" ");
    if (epos == std::string::npos) epos = word.find(";");
    word = word.substr(0, epos);
    if (! s2ul(word, &nchr)) {
        std::cout << "ERROR: Number of characters is not an integer"<< std::endl;
        exit(1);
    }

    // Set up split matrix
    for (int i = 0; i < nchr; i++) splits_.push_back("");

    // Roll forward until you find matrix
    found = false;
    while (std::getline(is, line)) {
        for (auto & c: line) c = std::tolower(c);
        sep = line.find("matrix");
        if (sep != std::string::npos) {
            found = true;
            break;
        }
    }
    if (!found) {
        std::cout << "ERROR: Unable to find matrix block in nexus file"<< std::endl;
        exit(1);
    }

    // Read data breaking when you get ;
    while (std::getline(is, line)) {
        sep = line.find(";");
        if (sep != std::string::npos) break;

        // Parse label
        sep = line.find(" ");
        if (sep == std::string::npos) {
            std::cout << "ERROR: Unable to read as phylip file"<< std::endl;
            exit(1);
        }
        word = line.substr(0, sep);
        labels_.push_back(word);
        //std::cout << "Added label " << word << std::endl;

        // Roll forward to sequence data
        for (pos = sep; pos < line.size(); pos++) {
            if (line[pos] != ' ') break;
        }
        sep = pos;

        if (line.size() - sep != (std::size_t) nchr) {
            std::cout << "ERROR: Wrong sequence length for " << word << std::endl;
            exit(1);
        }

        for (int i = 0; i < nchr; i++) {
            splits_[i].push_back(line[i + sep]);
        }
    }

    if (labels_.size() != (std::size_t) ntax) {
        std::cout << "ERROR: Read incorrect number of sequences" << std::endl;
        exit(1);
    }
}

void CharMat::read_phylip(std::istream &is) {
    // TODO: Pull in some error handling here:
    // https://github.com/ekmolloy/phylotools/blob/master/src/binary_character_matrix.cpp

    std::cout << "Reading input as phylip file" << std::endl;

    std::string line, word;
    std::size_t sep, pos;
    unsigned long ntax, nchr;

    // Read first line
    std::getline(is, line);

    sep = line.find(" ");
    if (sep == std::string::npos) {
         std::cout << "ERROR: Unable to parse first line"<< std::endl;
         exit(1);
    }
    
    // Parse into number of taxa
    word = line.substr(0, sep);
    if (! s2ul(word, &ntax)) {
        std::cout << "ERROR: Number of taxa is not an integer"<< std::endl;
        exit(1);
    }
    //std::cout << "Number of tax " << ntax << std::endl;

    // and number of characters
    word = line.substr(sep+1, std::string::npos);
    if (! s2ul(word, &nchr)) {
        std::cout << "ERROR: Number of characters is not an integer"<< std::endl;
        exit(1);
    }
    //std::cout << "Number of characters " << nchr << std::endl;

    // Set up split matrix
    for (int i = 0; i < nchr; i++) splits_.push_back("");

    while (std::getline(is, line)) {
        // Parse label
        sep = line.find(" ");
        if (sep == std::string::npos) {
            std::cout << "ERROR: Unable to read as phylip file"<< std::endl;
            exit(1);
        }
        word = line.substr(0, sep);
        labels_.push_back(word);
        //std::cout << "Added label " << word << std::endl;

        // Roll forward to sequence data
        for (pos = sep; pos < line.size(); pos++) {
            if (line[pos] != ' ') break;
        }
        sep = pos;

        if (line.size() - sep != (std::size_t) nchr) {
            std::cout << "ERROR: Wrong sequence length for " << word << std::endl;
            exit(1);
        }

        for (int i = 0; i < nchr; i++) {
            //std::cout << line[sep + i] << " x ";
            splits_[i].push_back(line[i + sep]);
        }
        //std::cout << std::endl;
    }

    if (labels_.size() != (std::size_t) ntax) {
        std::cout << "ERROR: Read incorrect number of sequences" << std::endl;
        exit(1);
    }
}
