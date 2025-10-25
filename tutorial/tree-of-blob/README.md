TOB-QMC Tutorial
==========================

1. Clone the TREE-QMC Github repository and build TREE-QMC following the instructions in the [README](../../README.md).


2. Go to tutorial directory.
```
cd tutorial/tree-of-blob
```

3. Run TOB-QMC on the [example input data](avian_uce_trees_3679.tre). The example data file contains the best maximum likelihood (ML) trees estimated for 3,679 UCEs with RAxML (note that branch support was estimated with rapid bootstrapping). This file comes from the [Avian Phylogenomics Project](https://doi.org/10.1186/s13742-014-0038-1).

*Below, we describe three different ways of running TREE-QMC given gene trees.*

Weighting quartets by branch support and branch length in the gene trees (hybrid mode)
---
Weighting quartets by branch support and branch length (hybrid mode) was originally proposed in [Weighted ASTRAL](https://doi.org/10.1093/molbev/msac215), referred to as **wASTRAL-h**. This approach can also be used with QTREE-QMC by using the `--hybrid` flag, along with one of the preset options to specify the type of branch support in the input. Try using the following command:
```
../../tree-qmc \
    --hybrid \
    --bootstrap \
    --root STRCA,TINMA \
    -i avian_uce_trees_3679.tre \
    -o treeqmc-hybrid.tre
```
which specifies the input gene trees have bootstrap support values and the species tree should be rooted at the clade containing `TINMA,STRCA` if possible. 
The `--support` option can be specified to annotate the branches of the species tree with quartet support. Alternatively, quartet support can be computed for a fixed species tree with the command:
```
../../tree-qmc \
    --hybrid \
    --bootstrap \
    --root STRCA,TINMA \
    --supportonly treeqmc-hybrid.tre \
    -i avian_uce_trees_3679.tre \
    -o annotated-treeqmc-hybrid.tre
```
Tip: Add `--writetable annotated-table.txt` to the command above if you would like the branch information and quartete support values written to a table.

Contracting low support branches in gene trees (instead of quartet weighting)
---
Contracting branches with low support is recommended when running [ASTRAL-III](https://doi.org/10.1186/s12859-018-2129-y). This approach can also be used with TREE-QMC, either by contracting branches in advance of running TREE-QMC *or* by using the `--contract <threshold>` flag. The command
```
../../tree-qmc \
    --bootstrap \
    --contract 0.10 \
    --root STRCA,TINMA \
    -i avian_uce_trees_3679.tre
```
specifies that branches in the input gene trees with bootstrap support less than `10` should be contracted (note: the contraction threshold is applied after mapping support values to the interval from 0 to 1).

Using the fast algorithm (no weighting or polytomies)
---
To run the fast TREE-QMC algorithm (no weighting or polytomies), use the command:
```
../../tree-qmc \
    --fast \
    --root STRCA,TINMA \
    -i avian_uce_trees_3679.tre
```
**IMPORTANT:** the above command does not allow quartet weighting and will cause any polytomies in the input to be randomly refined.
