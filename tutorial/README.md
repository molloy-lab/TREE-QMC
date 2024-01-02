Weighted TREE-QMC Tutorial
==========================

1. Clone the repository.
```
git clone https://github.com/molloy-lab/weighted-TREE-QMC.git
```

2. Build wTREE-QMC and go to tutorial directory.
```
cd weighted-TREE-QMC
cd MQLib
make
cd ..
g++ -std=c++11 -O2 -I MQLib/include -o wTREE-QMC src/*.cpp MQLib/bin/MQLib.a
```

3. Go to tutorial directory.
```
cd tutorial
```

3. Run TREE-QMC on the [example input data](avian_uce_trees_3679.tre). The example data file contains the best maximum likelihood (ML) trees estimated for 3,679 UCEs with RAxML (note that branch support was estimated with rapid bootstrapping). This file comes from the [Avian Phylogenomics Project](https://doi.org/10.1186/s13742-014-0038-1).

4. The output species tree is UNROOTED. **TODO** The root of the tree should be placed at `(TINMA,STRCA)`; see the provided [name map](avian_name_map.txt) for more information.

5. **TODO:** Estimate species tree branch support.

*Below, we describe four different ways of running steps 3 and 4.*

Weighting quartet by branch support only
---
Weighting quartets by branch support was originally proposed in [Weighted ASTRAL](https://doi.org/10.1093/molbev/msac215), referred to as **wASTRAL-s**. This approach can also be used with wQTREE-QMC (referred to as **wQTREE-QMC-s**) by using the `-w 1` flag, along with the `-r` flag to specify the minimum and maximum values that branch support can take on. For rapid bootstrapping with RAxML, the range should be set with `-r 0 100`.
```
../wTREE-QMC \
	-w 1 \
	-r 0 100 \
	-i avian_uce_trees_3679.tre \
	-o wtreeqmc-s.tre
```
**COMMON SETTINGS**: 
* For bootstrap support, set `-r 0 100`.
* For likelihood or probability support (e.g., sh), set `-r 0 1`.
* For local Bayesian support (e.g., abayes), set `-r 0.333 1`.

Weighting quartets by branch support and branch length (hybrid)
---
Weighting quartets by branch support and branch length (hybrid mode) was originally proposed in [Weighted ASTRAL](https://doi.org/10.1093/molbev/msac215), referred to as **wASTRAL-h**. This approach can also be used with wQTREE-QMC (referred to as **wQTREE-QMC-h**) by using the `-w 2` flag, along with the `-r` flag to specify the minimum and maximum values that branch support can take on.
```
../wTREE-QMC \
	-w 2 \
	-r 0 100 \
	-i avian_uce_trees_3679.tre \
	-o wtreeqmc-h.tre
```

Contracting low support branches (instead of quartet weighting)
---
Contracting branches with low support is recommended when running [ASTRAL-III](https://doi.org/10.1186/s12859-018-2129-y). This approach can also be used with wTREE-QMC, either by contracting branches in advance of running TREE-QMC and then using the `-w 0` flag *or* by using the `-c` flag, along with the support threshold. In the command below, all branches with bootstrap support less than `10` will be contracted.
```
../wTREE-QMC \
	-c 10 \
	-i avian_uce_trees_3679.tre \
	-o wtreeqmc-c10.tre
```

Weighting quartets by branch length only (not recommended)
---
Weighting quartets by branch length was originally proposed in [Weighted ASTRAL](https://doi.org/10.1093/molbev/msac215), referred to as **wASTRAL-bl**. This approach can also be used with wQTREE-QMC (referred to as **wQTREE-QMC-bl**) by using the `-w 3` flag. This option is **not recommended** but is provided for completeness.
```
../wTREE-QMC \
	-w 3 \
	-r 0 100 \
	-i avian_uce_trees_3679.tre \
	-o wtreeqmc-bl.tre
```

No quartet weighting
---
If you do not want to use quartet weighting, then you can use the faster TREE-QMC algorithm by specifying the `-w 4` command.
```
../wTREE-QMC \
	-w 4
	-i avian_uce_trees_3679.tre \
	-o wtreeqmc.tre
```
**Importantly**, if there are polytomies in the input, the above command will cause them to be randomly refined; random refinements can be avoided by using the slower `-w 0` option.
