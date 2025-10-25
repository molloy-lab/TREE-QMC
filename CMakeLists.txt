cmake_minimum_required(VERSION 3.18)
project(tree_qmc LANGUAGES CXX)

# ---------- compiler & build ----------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()


if(CMAKE_BUILD_TYPE STREQUAL "Release")
  add_compile_options(-O2 -DNDEBUG -fno-plt)
  add_link_options(-Wl,-O1,--as-needed)
endif()

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_ok OUTPUT ipo_msg)
if(ipo_ok)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
endif()

# ---------- self-contained R discovery (under external/R) ----------
set(R_ROOT "${CMAKE_SOURCE_DIR}/external/R" CACHE PATH "Self-contained R prefix")

set(_CANDIDATE_R_HOMES
  "${R_ROOT}/lib64/R"
  "${R_ROOT}/lib/R"
)
set(R_HOME "")
foreach(_p IN LISTS _CANDIDATE_R_HOMES)
  if(EXISTS "${_p}/include/R.h")
    set(R_HOME "${_p}")
    break()
  endif()
endforeach()
if(NOT R_HOME)
  message(FATAL_ERROR "Could not find R_HOME under ${R_ROOT} (expected lib64/R or lib/R with include/R.h)")
endif()

set(R_BIN        "${R_ROOT}/bin/R")
set(R_EXECUTABLE "${R_ROOT}/bin/Rscript")
if(NOT EXISTS "${R_BIN}")
  message(FATAL_ERROR "Missing R binary at ${R_BIN}")
endif()
if(NOT EXISTS "${R_EXECUTABLE}")
  message(FATAL_ERROR "Missing Rscript binary at ${R_EXECUTABLE}")
endif()

set(R_LIB "${R_HOME}/library" CACHE PATH "R site library (Rcpp, RInside, MSCquartets)")
if(NOT EXISTS "${R_HOME}/include/R.h")
  message(FATAL_ERROR "R headers not found at ${R_HOME}/include")
endif()
if(NOT EXISTS "${R_LIB}/RInside/include/RInside.h")
  message(FATAL_ERROR "RInside headers not found under ${R_LIB}/RInside/include")
endif()
if(NOT EXISTS "${R_LIB}/Rcpp/include/Rcpp.h")
  message(FATAL_ERROR "Rcpp headers not found under ${R_LIB}/Rcpp/include")
endif()
if(NOT EXISTS "${R_LIB}/MSCquartets")
  message(FATAL_ERROR "MSCquartets is not installed in ${R_LIB}. See bootstrap script.")
endif()

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
  "${R_HOME}/include"
  "${R_LIB}/Rcpp/include"
  "${R_LIB}/RInside/include"
)

# ---------- compile defs ----------
set(R_LIBDIRS "${R_LIB}")
target_compile_definitions(tree-qmc PRIVATE
  R_LIBDIRS="${R_LIBDIRS}"
  VERSION="42"
)

# ---------- link + rpath ----------
target_link_directories(tree-qmc PRIVATE
  "${R_HOME}/lib"               # libR.so
  "${R_LIB}/RInside/lib"        # libRInside.so
)
target_link_libraries(tree-qmc PRIVATE
  MQLib
  RInside
  R Rblas Rlapack
  m
)
set_target_properties(tree-qmc PROPERTIES
  BUILD_RPATH   "${R_HOME}/lib;${R_LIB}/RInside/lib"
  INSTALL_RPATH "${R_HOME}/lib;${R_LIB}/RInside/lib"
)

# Optional: place the binary elsewhere
# set_target_properties(tree-qmc PROPERTIES
#   RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/.."
# )
