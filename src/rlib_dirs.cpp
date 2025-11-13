// rlib_dirs.cpp

#include <string>
#include <iostream>
#include "rlib_dirs.hpp"

#ifndef R_LIBDIRS
#define R_LIBDIRS ""
#endif

#if TREE_QMC_WITH_R

// Add compile-time R library paths and load required R packages.
void add_r_libpaths_and_load(RInside& R) {
    const std::string libdirs = R_LIBDIRS;
    R[".tqmc_libdirs"] = libdirs;


    R.parseEvalQ(
        "suppressMessages({\n"
        "  if (nzchar(.tqmc_libdirs)) {\n"
        "    p <- strsplit(.tqmc_libdirs, .Platform$path.sep, fixed = TRUE)[[1]]\n"
        "    .libPaths(unique(c(p, .libPaths())))\n"
        "  }\n"
        "  if (!requireNamespace('MSCquartets', quietly = TRUE)) {\n"
        "    stop('MSCquartets not found. Current .libPaths() = ', paste(.libPaths(), collapse = ', '))\n"
        "  }\n"
        "  library(MSCquartets)\n"
        "})"
    );
}

#else

void add_r_libpaths_and_load(RInside& R) {
    // Do nothing if R is not enabled
}
#endif 