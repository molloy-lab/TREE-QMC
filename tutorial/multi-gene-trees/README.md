TREE-QMC-GDL Tutorial
==========================

1. Clone the TREE-QMC Github repository and build TREE-QMC following the instructions in the [README](../../README.md).

2. Go to tutorial directory.
```
cd tutorial/multi-gene-trees
```

3. Run TREE-QMC on the [example input data](multi_true_250_g_trees.trees).

*Below, we describe different ways of running TREE-QMC-GDL.*

Unrooted multi-labeled gene trees
---
TREE-QMC-GDL takes unrooted multi-labeled gene trees as input and performs rooting and tagging using the algorithm first introduced in [ASTRAL-Pro](https://doi.org/10.1093/molbev/msaa139). TREE-QMC-GDL uses the '--gdl' flag to be run using the duplication quartet exclusion principle. Normalization is not recomended in this mode so the '--norm_atax 0' option is required. Running with quartet weighting is supported but not tested, so we do not recommend running on weighted gene trees.
```
../../tree-qmc \
  --gdl \
  --norm-atax 0 \
  -i multi_true_250_g_trees.trees \
  -o treeqmc-unrooted.tre
```
the estimated species tree can be compared to the [true species tree](s_tree.trees).

Rooted and tagged multi-labeled gene trees
---
TREE-QMC-GDL can take rooted and tagged gene trees as input, allowing the user to use their preferred rooting and tagging algorithm. This approach can be used by using the '--tagged' flag in combination with the '--gdl' flag.
```
../../tree-qmc \
  --gdl \
  ---tagged \
  --norm-atax 0 \
  -i apro250_tagged_multi.trees \
  -o treeqmc-rooted-tagged.tre
```
