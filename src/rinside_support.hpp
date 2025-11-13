// rinside_support.hpp
#pragma once

#ifndef TREE_QMC_WITH_R
#define TREE_QMC_WITH_R 0
#endif

#if TREE_QMC_WITH_R
  #include <RInside.h>
#else
  // Dummy RInside when R is disabled.
  class RInside {};
#endif
