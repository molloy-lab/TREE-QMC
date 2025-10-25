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
    if (!libdirs.empty()) {
        // allow multiple colon-separated paths
        R.parseEvalQ(
            "suppressMessages({"
            "  p <- strsplit(\"" + libdirs + "\", \":\")[[1]];"
            "  .libPaths(c(p, .libPaths()))"
            "})"
        );
    }

    // Load MSCquartets and TINNIK
    R.parseEvalQ(
        "suppressMessages({"
        "  tryCatch(library(MSCquartets), error=function(e){stop('MSCquartets not found: ', e$message)})"
        "})"
    );
}
