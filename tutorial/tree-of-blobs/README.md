TOB-QMC Tutorial
==========================

1. Clone the TREE-QMC Github repository and build TREE-QMC following the instructions in the [README](../../README.md).


2. Go to tutorial directory.
```
cd tutorial/tree-of-blobs
```

3. Run TOB-QMC on the [example input data](gene_trees.tre). The example data file contains 991 the best maximum likelihood (ML) trees estimated with iqtree3 (note that branch support was estimated with rapid bootstrapping).

To infer a tree of blobs, our approach requires frstly build a base tree(a refinement of the tree of blobs)

*Below, we describe three different ways of running TOB-QMC given gene trees to get a base tree.*

Build a base tree by hill-climbing heuristic
---
This approach can be used with TOB-QMC by using the `--store_pvalue --iter_limit_blob` flags, along with one of the preset options to specify the type of branch support in the input. Try using the following command:
```
../../build/tree-qmc \
    --iter_limit_blob 5000 \
    --store_pvalue \
    -i gene_trees.tre \
    -o base_tree.tre
```
which specifies the maximum number of iterations applied in hill-climbing heurstic for each branch. By default, we recommand set it to be $2n^2$, where $n$ is the number of taxa. In this example it is $5000$. 
Tip: Add `--override` to the command above if you would like to try different seeting for the maximum number of iterations and want to overwrite the output base tree file.
Tips: Add `--override` to the above command above if you would like to explore different maximum number of iterations and overwrite the output file.

Build a base tree by exhaustive all quartets
---
This approach can be used with TREE-QMC, by simply set `--iter_limit_blob 0` . The command
```
../../build/tree-qmc \
    --iter_limit_blob 0 \
    --store_pvalue \
    -i gene_trees.tre \
    -o base_tree_exhaustive.tre
```

Using the 3-fix-1-alter algorithm
---
This approach will use 3-fix-1-alter to build the base tree. This approach is statistical consistent and only sampling $\Theta(n)$ quartets for each branch. To run the 3-fix-1-alter TOB-QMC algorithm, use the command:
```
../../build/tree-qmc \
    --3f1a \
    --store_pvalue \
    -i gene_trees.tre \
    -o base_tree_3f1a.tre
```
**IMPORTANT:** We only recommand use this algorithm all quartet concordance factors are close to its expectation.

Once get a base tree, TOB-QMC can contract branches in the base tree to infer a tree of blobs, use the command:
```
../../build/tree-qmc \
    --blob \
    --alpha 1e-7 \
    --beta 0.95 \
    --load_pvalue \
    -i base_tree.tre \
    -o tob_default.tre \
```
which species the hyperparameter $\alpha$ and $\beta$ used in hypothesis testing be $10^{-7}$ and $0.95$, respectively. 
Tips: Add `--override` to the above command above if you would like to explore different setting of hyperparameters and overwrite the output file.