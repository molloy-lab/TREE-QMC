#if ENABLE_TOB

// rlib_dirs.cpp
#include <RInside.h>
#include <string>
#include <iostream>

#ifndef R_LIBDIRS
#define R_LIBDIRS ""
#endif

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
#endif  // ENABLE_TOB