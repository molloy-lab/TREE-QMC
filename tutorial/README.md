Weighted TREE-QMC Tutorial
==========================

1. Clone the repository.
```
git clone https://github.com/molloy-lab/weighted-TREE-QMC.git
```

2. Build wTREE-QMC and go to tutorial directory.
```
git clone https://github.com/molloy-lab/weighted-TREE-QMC.git
cd weighted-TREE-QMC/external/MQLib
make
cd ../..
g++ -std=c++11 -O2 \
    -I external/MQLib/include -I external/toms743 \
    -o treeqmc \
    src/*.cpp external/toms743/toms743.cpp \
    external/MQLib/bin/MQLib.a -lm 
```

3. Go to tutorial directory.
```
cd tutorial
```

4. Run TREE-QMC on the [example input data](avian_uce_trees_3679.tre). The example data file contains the best maximum likelihood (ML) trees estimated for 3,679 UCEs with RAxML (note that branch support was estimated with rapid bootstrapping). This file comes from the [Avian Phylogenomics Project](https://doi.org/10.1186/s13742-014-0038-1).

*Below, we describe four different ways of running steps 3 and 4.*

Weighting quartets by branch support and branch length (hybrid)
---
Weighting quartets by branch support and branch length (hybrid mode) was originally proposed in [Weighted ASTRAL](https://doi.org/10.1093/molbev/msac215), referred to as **wASTRAL-h**. This approach can also be used with QTREE-QMC by using the `--hybrid` flag, along with one of the preset options to specify the type of branch support in the input. Try using the following command:
```
../treeqmc \
	--hybrid \
	--bootstrap \
	--root STRCA,TINMA \
	-i avian_uce_trees_3679.tre \
	-o treeqmc-hybrid.tre
```
which specifies the input gene trees have bootstrap support values and the species tree should be rooted at the clade containing `TINMA,STRCA` if possible. 
The `--support` option can be specified to annotate the branches of the species tree with quartet support. Alternatively, quartet support can be computed for a fixed species tree with the command:
```
../treeqmc \
	--hybrid \
	--bootstrap \
	--supportonly treeqmc-hybrid.tre \
	-i avian_uce_trees_3679.tre \
	-o annotated-treeqmc-hybrid.tre
```
Lastly, the `--writetable <output file name>` option can be included in the command above to additionally write the branch support information and support values written to a table.

Contracting low support branches (instead of quartet weighting)
---
Contracting branches with low support is recommended when running [ASTRAL-III](https://doi.org/10.1186/s12859-018-2129-y). This approach can also be used with TREE-QMC, either by contracting branches in advance of running TREE-QMC *or* by using the `--contract <threshold>` flag. Try using the command:
```
../treeqmc \
	--bootstrap \
	--contract 0.10 \
	--root STRCA,TINMA \
	-i avian_uce_trees_3679.tre
```
which specifies that branches in the input gene trees with bootstrap support less than `10` should be contracted (note: the contraction threshold is applied after mapping support values to the interval from 0 to 1).

Using the faster algorithm
---
The faster TREE-QMC algorithm (no weighting) can invoked with the command:
```
../treeqmc \
	--fast \
	-i avian_uce_trees_3679.tre
```
**Importantly**, the above command does not allow quartet weighting and will cause any polytomies in the input to be randomly refined.
