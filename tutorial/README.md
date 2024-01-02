Weighted TREE-QMC Tutorial
==========================

1. Clone the repository.
```
git clone https://github.com/molloy-lab/weighted-TREE-QMC.git
```

2. Build wTREE-QMC.
```
cd weighted-TREE-QMC
cd MQLib
make
cd ..
g++ -std=c++11 -O2 -I MQLib/include -o wTREE-QMC src/*.cpp MQLib/bin/MQLib.a
```

3. Run TREE-QMC on the [example input data](avian_uce_trees_3679.tre). The example data file contains the best maximum likelihood (ML) trees estimated for 3,679 UCEs with RAxML (note that branch support was estimated with rapid bootstrapping). This file comes from the [Avian Phylogenomics Project](https://doi.org/10.1186/s13742-014-0038-1).

4. Estimate species tree branch support using either ASTRAL or Weighted ASTRAL/ASTER. In the future, we plan to implement these calculations within TREE-QMC for ease of use.

IMPORTANT: The output species tree is UNROOTED. The root of the tree should be placed at `(TINMA,STRCA)`; see the provided [name map](avian_name_map.txt) for more information.

Below, we describe four different ways of running steps 3 and 4. 

