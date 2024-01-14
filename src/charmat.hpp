#ifndef CHARMAT_HPP
#define CHARMAT_HPP

#include "utility.hpp"

class CharMat {
	public:
        CharMat(std::string input_file);
        ~CharMat();
        std::size_t size();
        void pop_newick(std::string &newick);
	private:
        std::vector<std::string> splits_;  // actually characters...
        std::vector<std::string> labels_;
        void read_phylip(std::istream &is);
        void read_fasta(std::istream &is);
        void read_nexus(std::istream &is);
};

#endif
