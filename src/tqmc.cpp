#ifndef ENABLE_TOB
#define ENABLE_TOB FALSE
#endif

#include "instance.hpp"

std::ofstream subproblem_csv, quartets_txt, good_edges_txt, bad_edges_txt;
std::string verbose = "0";
unsigned long long count[8];

int main(int argc, char **argv) {
    auto start = std::chrono::high_resolution_clock::now();

    Instance instance(argc, argv);
    long long time = instance.solve();

    if (instance.get_solution() != NULL)
        instance.output_solution();

    auto end = std::chrono::high_resolution_clock::now();

    //const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    //std::cout << "Execution time: " << ms.count() << "ms" << std::endl;

    const auto secs = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "Execution time: " << secs.count() << "secs" << std::endl;

    return 0;
}
