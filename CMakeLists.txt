cmake_minimum_required(VERSION 3.18)
project(tree_qmc LANGUAGES CXX)

# ---------- compiler & build ----------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(APPLE)
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
	add_compile_options(-O2 -DNDEBUG)
        add_link_options(-Wl,-dead_strip)
    endif()
else()
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
	add_compile_options(-O2 -DNDEBUG -fno-plt)
        add_link_options(-Wl,-O1,--as-needed)
    endif()
endif()

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_ok OUTPUT ipo_msg)
if(ipo_ok)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
endif()

# ---------- system R discovery (replaces external/R) ----------
find_program(RSCRIPT_BIN NAMES Rscript
  HINTS /opt/homebrew/bin /usr/local/bin /usr/bin
  REQUIRED)
message(STATUS "Found Rscript at: ${RSCRIPT_BIN}")

# Prefer env overrides; otherwise query R
if(DEFINED ENV{R_INCLUDE_DIR})
  set(R_INCLUDE_DIR "$ENV{R_INCLUDE_DIR}")
else()
  execute_process(
    COMMAND "${RSCRIPT_BIN}" --vanilla -e "cat(normalizePath(R.home('include')))"
    OUTPUT_VARIABLE R_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
endif()

execute_process(
  COMMAND "${RSCRIPT_BIN}" --vanilla -e "cat(normalizePath(R.home()))"
  OUTPUT_VARIABLE R_HOME
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(DEFINED ENV{R_LIB_DIR})
  set(R_LIB_DIR "$ENV{R_LIB_DIR}")
else()
  execute_process(
    COMMAND "${RSCRIPT_BIN}" --vanilla -e "cat(normalizePath(R.home('lib')))"
    OUTPUT_VARIABLE R_LIB_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
endif()

# Validate headers
if(NOT R_INCLUDE_DIR OR NOT EXISTS "${R_INCLUDE_DIR}/R.h")
  message(FATAL_ERROR
    "Could not locate R headers. Looked at R_INCLUDE_DIR='${R_INCLUDE_DIR}'.\n"
    "Tip: export R_INCLUDE_DIR=/usr/share/R/include")
endif()

message(STATUS "R_HOME       = ${R_HOME}")
message(STATUS "R_INCLUDE_DIR= ${R_INCLUDE_DIR}")
message(STATUS "R_LIB_DIR    = ${R_LIB_DIR}")

# RInside + Rcpp dirs
execute_process(
  COMMAND "${RSCRIPT_BIN}" --vanilla -e "p<-system.file('include',package='RInside');cat(if(nzchar(p)) normalizePath(p) else '')"
  OUTPUT_VARIABLE RINSIDE_INC
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT EXISTS "${RINSIDE_INC}/RInside.h")
  message(FATAL_ERROR "RInside headers not found for this R. In R: install.packages('RInside')")
endif()

execute_process(
  COMMAND "${RSCRIPT_BIN}" --vanilla -e "p<-system.file('lib',package='RInside');cat(if(nzchar(p)) normalizePath(p) else '')"
  OUTPUT_VARIABLE RINSIDE_LIB_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT EXISTS "${RINSIDE_LIB_DIR}")
  message(FATAL_ERROR "RInside lib dir not found. In R: install.packages('RInside')")
endif()

execute_process(
  COMMAND "${RSCRIPT_BIN}" --vanilla -e "p<-system.file('include',package='Rcpp');cat(if(nzchar(p)) normalizePath(p) else '')"
  OUTPUT_VARIABLE RCPP_INC
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT EXISTS "${RCPP_INC}/Rcpp.h")
  message(FATAL_ERROR "Rcpp headers not found. In R: install.packages('Rcpp')")
endif()

# ensure MSCquartets installed
execute_process(
  COMMAND "${RSCRIPT_BIN}" --vanilla -e "cat(if(requireNamespace('MSCquartets',quietly=TRUE)) 'OK' else '')"
  OUTPUT_VARIABLE MSCQ_OK
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT MSCQ_OK STREQUAL "OK")
  message(FATAL_ERROR "R package 'MSCquartets' not found. In R: install.packages('MSCquartets')")
endif()

# --- Gather R's library trees for runtime (used by your RInside code) ---
execute_process(
  COMMAND "${RSCRIPT_BIN}" --vanilla -e "cat(paste(normalizePath(.libPaths()), collapse=':'))"
  OUTPUT_VARIABLE R_LIB_TREES
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
string(REPLACE "\"" "\\\"" R_LIB_TREES "${R_LIB_TREES}")
message(STATUS "R library trees: ${R_LIB_TREES}")

# ---------- external MQLib (built via its Makefile) ----------
set(MQLIB_DIR "${CMAKE_SOURCE_DIR}/external/MQLib")
set(MQLIB_A   "${MQLIB_DIR}/bin/MQLib.a")
if(NOT EXISTS "${MQLIB_DIR}/Makefile")
  message(FATAL_ERROR "Missing ${MQLIB_DIR}/Makefile")
endif()

add_custom_command(
  OUTPUT  ${MQLIB_A}
  COMMAND ${CMAKE_COMMAND} -E env CXXFLAGS="-O2 -std=c++17"
          ${CMAKE_MAKE_PROGRAM} -C ${MQLIB_DIR}
  WORKING_DIRECTORY ${MQLIB_DIR}
  COMMENT "Building MQLib with make (-O2 -std=c++17)"
  VERBATIM
)
add_custom_target(build_mqlib ALL DEPENDS ${MQLIB_A})

add_library(MQLib STATIC IMPORTED GLOBAL)
set_target_properties(MQLib PROPERTIES IMPORTED_LOCATION ${MQLIB_A})
add_dependencies(MQLib build_mqlib)

# ---------- sources ----------
file(GLOB TREE_QMC_SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")
if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/toms743/toms743.cpp")
  message(FATAL_ERROR "Missing toms743.cpp at external/toms743/toms743.cpp")
endif()
list(APPEND TREE_QMC_SOURCES "${CMAKE_SOURCE_DIR}/external/toms743/toms743.cpp")
add_executable(tree-qmc ${TREE_QMC_SOURCES})

# ---------- includes ----------
target_include_directories(tree-qmc
  SYSTEM PRIVATE
  "${CMAKE_SOURCE_DIR}/external/MQLib/include"
)
target_include_directories(tree-qmc PRIVATE
  "${CMAKE_SOURCE_DIR}/src"
  "${CMAKE_SOURCE_DIR}/external/toms743"
  "${R_INCLUDE_DIR}"
  "${RCPP_INC}"
  "${RINSIDE_INC}"
)

# ---------- compile defs ----------
file(READ "${CMAKE_SOURCE_DIR}/version.txt" VERSION_TXT)
string(STRIP "${VERSION_TXT}" VERSION_TXT)
target_compile_definitions(tree-qmc PRIVATE VERSION="${VERSION_TXT}" -DENABLE_TOB)
target_compile_definitions(tree-qmc PRIVATE R_LIBDIRS="${R_LIB_TREES}" -DENABLE_TOB)

# ---------- link + rpath ----------
# Link search paths (only add if they exist)
if(EXISTS "${R_LIB_DIR}")
  target_link_directories(tree-qmc PRIVATE "${R_LIB_DIR}")
endif()
if(EXISTS "${RINSIDE_LIB_DIR}")
  target_link_directories(tree-qmc PRIVATE "${RINSIDE_LIB_DIR}")
endif()

if(APPLE)
  # Resolve actual dylib files
  set(RINSIDE_DYLIB "${RINSIDE_LIB_DIR}/libRInside.dylib")
  if(NOT EXISTS "${RINSIDE_DYLIB}")
    message(FATAL_ERROR "Expected RInside dylib not found at: ${RINSIDE_DYLIB}")
  endif()

  # libR.dylib candidates
  set(R_DYLIB_CANDIDATES
    "${R_LIB_DIR}/libR.dylib"
    "${R_HOME}/lib/libR.dylib"
    "${R_HOME}/lib/R/lib/libR.dylib"
  )
  set(R_DYLIB "")
  foreach(_c IN LISTS R_DYLIB_CANDIDATES)
    if(EXISTS "${_c}")
      set(R_DYLIB "${_c}")
      break()
    endif()
  endforeach()
  if(NOT R_DYLIB)
    message(FATAL_ERROR "Could not find libR.dylib under expected locations.")
  endif()

  target_link_libraries(tree-qmc PRIVATE
    MQLib
    "${RINSIDE_DYLIB}"
    "${R_DYLIB}"
  )
else()
  target_link_libraries(tree-qmc PRIVATE
    MQLib
    RInside
    R
    m
  )
  # Keep the linker from dropping libR/libRInside under --as-needed
  target_link_options(tree-qmc PRIVATE -Wl,--no-as-needed)
endif()

# RPATHs
if(APPLE)
  set(_RPATHS "${RINSIDE_LIB_DIR};${R_LIB_DIR};@loader_path;@loader_path/../lib;@loader_path/../R/lib")
  set_target_properties(tree-qmc PROPERTIES
    BUILD_RPATH   "${_RPATHS}"
    INSTALL_RPATH "${_RPATHS}"
  )
else()
  # Make runtime find libRInside.so and libR.so
  set(_LINUX_RPATHS "${RINSIDE_LIB_DIR};${R_LIB_DIR}")
  set_target_properties(tree-qmc PROPERTIES
    BUILD_RPATH   "${_LINUX_RPATHS}"
    INSTALL_RPATH "${_LINUX_RPATHS}"
  )
endif()

# macOS post-link: only rewrite install names (no rpath additions here)
if(APPLE)
  add_custom_command(TARGET tree-qmc POST_BUILD
    COMMAND /usr/bin/install_name_tool -change libRInside.dylib "${RINSIDE_DYLIB}" $<TARGET_FILE:tree-qmc>
    COMMAND /usr/bin/install_name_tool -change libR.dylib      "${R_DYLIB}"       $<TARGET_FILE:tree-qmc>
    VERBATIM
  )
endif()

# Optional: place the binary elsewhere
# set_target_properties(tree-qmc PROPERTIES
#   RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/.."
# )
