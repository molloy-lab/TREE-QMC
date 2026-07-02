#include "network.hpp"

NetworkNode::NetworkNode(std::string label) {
    this->label = label;
    this->length = 0;
    this->support = 1;
    this->lambda = -1;
    this->hybrid = this->parent = NULL;
    count ++;
}

NetworkNode::~NetworkNode() {
    count --;
    for (auto child : children) 
        delete child;
}

std::size_t NetworkNode::count = 0;

std::size_t NetworkNode::get_count() {
    return count;
}
