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
    std::size_t ntax, ntax_done, ntax_found;
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
            itr->second.first += "," + labels_[j];
            itr->second.second++;
            //std::cout << itr->first << ": " << itr->second.first << " (" << itr->second.second << " taxa)" << std::endl;
        }
    }

    // Convert to newick string
    newick = "(";
    ntax_done = 0;

    // Remove missing taxa
    auto itr = mymap.find('?');
    if (itr != mymap.end()) ntax_done += itr->second.second;

    // Process ancestral state
    itr = mymap.find('0');
    if (itr != mymap.end()) {
        newick += itr->second.first;
        ntax_done += itr->second.second;
        if (ntax_done != ntax) newick += ",";
    }

    for (itr = mymap.begin(); itr != mymap.end(); ++itr) {
        if (itr->first != '?' && itr->first != '0') {
            ntax_found = itr->second.second;

            if (ntax_found == 1)
                newick += itr->second.first;
            else
                newick += "(" + itr->second.first + ")";

            ntax_done += ntax_found;
            if (ntax_done != ntax) newick += ",";
        }
    }

    newick += ");";

    //std::cout << split << std::endl;
    //std::cout << newick << std::endl;

    // Pop split
    splits_.pop_back();
}

index_t CharMat::size() {
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

    // Count characters and set up split matrix
    std::getline(is, line);
    if (line[0] == '>') {
        std::cout << "ERROR: Unable to read as fasta file" << std::endl;
        exit(1);
    }
    nchr = line.size();
    for (int i = 0; i < nchr; i++) {
        splits_.push_back("");
        splits_[i].push_back(line[i]);
    }
    nextlab = true;

    // Read data
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
}

void CharMat::read_nexus(std::istream &is) {
    std::cout << "Reading input as nexus file" << std::endl;
    exit(1);
}

void CharMat::read_phylip(std::istream &is) {
    std::cout << "Reading input as phylip file" << std::endl;
    exit(1);
}


