#ifndef DICT_HPP
#define DICT_HPP

#include "utility.hpp"

class Dict {
    public:
        index_t label2index(std::string label) {
            if (l2i.find(label) == l2i.end()) {
                l2i.insert(std::make_pair(label, i2l.size()));
                i2l.push_back(label);
            }
            return l2i[label];
        }

        std::string index2label(index_t index) {
            if (index < 0 || index >= singletons) 
                return "A" + std::to_string(index);
            return i2l[index];
        }

        index_t size() {
            return i2l.size();
        }

        std::string to_string() {
            std::string s = "";
            for (std::string label : i2l) 
                s += label + "->" + std::to_string(l2i[label]) + "\n";
            return s;
        }

        void update_singletons() {
            singletons = size();
        }

        index_t max_size() {
            return singletons * 2 - 3;
        }

    private:
        std::unordered_map<std::string, index_t> l2i;
        std::vector<std::string> i2l;
        index_t singletons;
};

#endif
