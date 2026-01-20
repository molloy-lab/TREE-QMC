#include "utility.hpp"

bool DEBUG_MODE;

weight_t **Matrix::new_mat(index_t size) {
    weight_t **m = new weight_t*[size];
    for (index_t i = 0; i < size; i ++) {
        m[i] = new weight_t[size];
        for (index_t j = 0; j < size; j ++) {
            m[i][j] = 0;
        }
    }
    return m;
}

void Matrix::delete_mat(weight_t **m, index_t size) {
    for (index_t i = 0; i < size; i ++) {
        delete[] m[i];
    }
    delete[] m;
}

std::string Matrix::display_mat(weight_t **m, index_t size) {
    std::stringstream ss;
    for (index_t i = 0; i < size; i ++) {
        for (index_t j = 0; j < size; j ++) {
            ss << std::setw(12) << std::setprecision(6) << (double)m[i][j];
        }
        ss << std::endl;
    }
    return ss.str();
}

weight_t Matrix::diff_mat(weight_t **m1, weight_t **m2, index_t size) {
    weight_t sum = 0;
    for (index_t i = 0; i < size; i ++) {
        for (index_t j = 0; j < size; j ++) {
            weight_t delta = m1[i][j] - m2[i][j];
            if (delta < 0) delta = - delta;
            sum += delta;
        }
    }
    return sum;
}

quartet_t join(index_t *indices) {
    index_t temp;
    if (indices[0] < indices[1]) {
        temp = indices[0]; 
        indices[0] = indices[1]; 
        indices[1] = temp;
    }
    if (indices[2] < indices[3]) {
        temp = indices[0]; 
        indices[0] = indices[1]; 
        indices[1] = temp;
    }
    if ((quartet_t) indices[0] * INDEX_WIDTH + indices[1] < (quartet_t) indices[2] * INDEX_WIDTH + indices[3]) {
        temp = indices[0]; 
        indices[0] = indices[2]; 
        indices[2] = temp;
        temp = indices[1]; 
        indices[1] = indices[3]; 
        indices[3] = temp;
    }
    quartet_t quartet = 0;
    for (index_t i = 0; i < 4; i ++) 
        quartet = quartet * INDEX_WIDTH + indices[i];
    return quartet;
}

index_t *split(quartet_t quartet) {
    index_t *indices = new index_t[4];
    for (index_t i = 3; i >= 0; i --) {
        indices[i] = quartet % INDEX_WIDTH;
        quartet /= INDEX_WIDTH;
    }
    return indices;
}

bool s2d(std::string s, weight_t *r) {
    char* end;
    weight_t tmp = strtod(s.c_str(), &end);
    if (end == s.c_str() || *end != '\0') return false;
    *r = tmp;
    return true;
}

bool s2ul(std::string s, unsigned long int *r) {
    char* end;
    unsigned long int tmp = strtoul(s.c_str(), &end, 10);
    if (end == s.c_str() || *end != '\0') return false;
    *r = tmp;
    return true;
}

weight_t* init(index_t size) {
    weight_t *a = new weight_t[size];
    for (index_t i = 0; i < size; i ++) 
        a[i] = std::nan("");
    return a;
}
