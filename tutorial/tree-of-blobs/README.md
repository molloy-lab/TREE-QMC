TOB-QMC Tutorial
=================

This tutorial shows how to run TREE-QMC to reconstruct a tree of blobs (TOB), which requires TREE-QMC to be built with R.

---

BUILD
-----

Requirements
* [cmake](https://cmake.org/download/)
* [R](https://www.r-project.org/)
* [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)
* [RInside](https://cran.r-project.org/web/packages/RInside/index.html) 
* [MSCquartets](https://cran.r-project.org/web/packages/MSCquartets/)

**Step 1.** If cmake isn't already installed on your system, install cmake version 3.18 or later by following the instructions on the [cmake official website](https://cmake.org/download/). If you are using a cluster, you may be able to load cmake as a module. Type `module avail cmake` to see the available packages and then load cmake, if available, with the `module load <package>` command.

**Step 2.** If R isn't already installed on your system, install R following the instructions on the [R official website](https://www.r-project.org/). If you are using a cluster, you may be able to load R as a module. Type `module avail R` to see the available packages and then load R, if available, with the `module load <package>` command. 

**Step 3.** Start the R console by typing `R` into the terminal as a commandline instruction and then install packages by typing
```
install.packages('Rcpp', repos='https://cloud.r-project.org/')
install.packages('RInside', repos='https://cloud.r-project.org/')
install.packages('MSCquartets', repos='https://cloud.r-project.org/')
```
If you are working on a cluster, you may get the following Warning message:
```
Warning in install.packages("Rcpp", repos = "https://cloud.r-project.org/") :
  'lib = "/opt/local/stow/R-4.5.1/lib64/R/library"' is not writable
Would you like to use a personal library instead? (yes/No/cancel) yes
Would you like to create a personal library
‘~/R/x86_64-pc-linux-gnu-library/4.5’
to install packages into? (yes/No/cancel)
```
Just type `yes` to install packages into a personal library.

**Step 4.** Build TREE-QMC.
```
git clone https://github.com/molloy-lab/TREE-QMC.git
cd TREE-QMC
mkdir -p build
cd build
cmake ..              
cmake --build . -j
```

**Step 5.** Add TREE-QMC build directory to your path.
```
TREEQMC_PATH=$(pwd)
export PATH="$TREEQMC_PATH:$PATH"
cd ..
```

TUTORIAL
-----
**Step 1.** Go to tutorial directory.
```
cd tutorial/tree-of-blobs
```

**Step 2.** Run TOB-QMC on the [example input data](nomiinae_gene_trees.txt). The example data file contains gene trees for the bee subfamily *Nomiinae* estimated for 852 UCEs.

To infer a tree of blobs, our approach requires that you first build a base tree with branches annotated by the minimum p-value found from applying hypothesis tests to 4-taxon subsets around each edge. After, edges are contracted based on the alpha and beta hyperparameter values.

*Below, we describe three different ways of running TOB-QMC to get a base tree.*

Build a base tree using the bipartition search heuristic.
---
This approach can be used with TOB-QMC by using the `--store_pvalue --iter_limit_blob` flags, along with one of the preset options to specify the type of branch support in the input. Try using the following command:
```
tree-qmc \
    --iter_limit_blob 2178 \
    --store_pvalue \
    --root Lasioglossum_albipes \
    -i nomiinae_gene_trees.tre \
    -o nomiinae_base_tree_default.tre
```
which specifies the maximum number of iterations for each branch. By default, we recommand set it to be $2n^2$, where $n$ is the number of taxa. In this example, it is $1922$ because there are $31$ taxa.

**Tip:** Add `--override` to the command above or change the output file name if you would like to try different settings for the maximum number of iterations and want to overwrite the output base tree file.

Build a base tree by exhaustive testing all 4-taxon subsets around each branch (bipartition)
---
This approach can be used with TREE-QMC, by simply set `--iter_limit_blob 0`. Try using the following command:
```
tree-qmc \
    --iter_limit_blob 0 \
    --store_pvalue \
    --root Lasioglossum_albipes \
    -i nomiinae_gene_trees.tre \
    -o nomiinae_base_tree_exhaustive.tre
```

Build a base tree by using the 3-fix-1-alter (3f1a) search heuristic.
---
The 3f1a search algorithm is statistically consistent and only samples $\Theta(n)$ quartets for each branch. To run the 3-fix-1-alter algorithm, use the command:
```
tree-qmc \
    --3f1a \
    --store_pvalue \
    --root Lasioglossum_albipes \
    -i nomiinae_gene_trees.tre \
    -o nomiinae_base_tree_3f1a.tre
```
**IMPORTANT:** We only recommand use this algorithm all quartet concordance factors are close to their expectation. In our experimental study, the default (bipartition search heuristic) outperforms the 3f1a algorithm.

Contract branches based on hyperparameter thresholds.
---
Once you have the base tree, TOB-QMC can contract branches based on the alpha and beta hyperparameters to infer a tree of blobs. Use the command:
```
tree-qmc \
    --blob \
    --alpha 1e-6 \
    --beta 0.95 \
    --load_pvalue \
    -i nomiinae_base_tree_default.tre \
    -o nomiinae_tob_default.tre
```
which sets the $\alpha$ and $\beta$ hyperparameters to $10^{-6}$ and $0.95$, respectively.
This contracts one branch for the genus Stictonomia based on the following tree-test:
```
QTT: 1.36381e-07 [69/80/21] minimizer: [Nomiapis_bispinosa/Stictonomia_sangaensis/Stictonomia_aliceae/Stictonomia_schubotzi]
```
This line in the logfile from constructing the base tree shows that the min p-value found around the branch is 1.36381e-07. This p-value occurs for quartet concordance factors (qCFs) equal to 69, 80, and 21. These qCFs are suggestive of non-tree-like evolution (otherwise we expect the two lower qCFs to be equal to each other). That being said, these four taxa appear in only 170 out of the 852 UCE trees so there is limited data for hypothesis testing.

**Tip:** Add `--override` to the command above or change the output file name if you would like to try different settings for the maximum number of iterations and want to overwrite the output base tree file.
