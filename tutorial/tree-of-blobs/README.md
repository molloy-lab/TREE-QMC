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
```

TUTORIAL
-----
In this tutorial, we will reconstruct a tree of blobs (TOB) for the bee subfamily *Nomiinae* (31 taxa) given [852 UCE gene trees](nomiinae_gene_trees.txt). To begin, go to tutorial directory.
```
cd ../tutorial/tree-of-blobs/nomiinae
```
To run TOB-QMC in default mode, type:
```
tree-qmc \
    --blob \
    -i nomiinae_gene_trees.tre \
    -o tob_default_nominaee.tre
```
*Tip:* Add `--override` to the command above or change the output file name.

The command above reconstructs a TOB in three steps:
1. build a base tree (also called a refinement tree)
2. search for network-like-signal around each edge in the base tree (i.e., search for mininum p-value by applying TINNiK's tree and star tests to 4-taxon subsets around each edge in the base tree)
3. contract edges in the base tree applying significance levels for hypothesis testing. 

However, the output TOB is binary (try plotting it in [IcyTree](https://icytree.org)), and you don't learn much about your data. Instead, it is helpful to run each step of TOB-QMC independently. Instead, it can be helpful to execute each step separately

Step 1: Build a base tree (aka refinement TOB) with TREE-QMC
---
```
../../../build/tree-qmc \
    --root Lasioglossum_albipes \
    -i nomiinae_gene_trees.tre \
    -o nomiinae_base_tree.tre
```

Step 2: Annotate branches with minimum p-value found from hypothesis testing
---
There are several different heuristics you can use to search for a minimum p-value. To execute the default search, type
```
tree-qmc \
    --blobsearchonly nomiinae_base_tree.tre \
    --store_pvalue \
    -i nomiinae_gene_trees.tre \
    -o nomiinae_psearch_default.tre
```
This command executes bipartition search, which samples 2 taxa on either side of a branch. The default iteration limit (i.e., maximum number of 4-taxon subsets to test) is two times the number of species squared; this can be changed by adding the flag `--iter_limit_blob <num>`. For large data sets, we recommend setting the iteration limit to one fourth the number of taxa squared. For small data sets, we recommend exhaustive search around the branch (bipartition). To execute exhaustive search, type
```
tree-qmc \
    --blobsearchonly nomiinae_base_tree.tre \
    --store_pvalue \
    --iter_limit_blob 0 \
    -i nomiinae_gene_trees.tre \
    -o nomiinae_psearch_exhaustive.tre \
    &> nomiinae_psearch_exhaustive.log
```
Then type
```
cat nomiinae_psearch_exhaustive.log
```
to look at the log file.

Step 3: Contract branches based on hyperparameter settings
---
A major challenge with hypothesis testing in this context is selecting the hyperparameters alpha and beta for the quartet tree test (QTT) and the quartet star test (QST), respectively. However, if you look at the log file, it appears that there is one branch with some signal of hybridization; specifically, line
```
QTT: 1.36381e-07 [69/80/21] minimizer: [Nomiapis_bispinosa/Stictonomia_sangaensis/Stictonomia_aliceae/Stictonomia_schubotzi]
```
indicates a min p-value of 1.36381e-07 for qCFs of 69, 80, and 21. These qCFs are suggestive of non-tree-like evolution because the two smaller values are not equal to each other. To contract the branch, type
```
tree-qmc \
    --alpha 1e-6 \
    --beta 0.95 \
    --load_pvalue \
    -i nomiinae_psearch_exhaustive.tre \
    -o nomiinae_psearch_exhaustive_a1e-6_b0.95_tob.tre
```
The final TOB has the genus Stictonomia as a blob. Lastly, we confirmed the blob signal by testing all other 4-taxon subsets around the branch; see [script](test_around_branch_with_TINNiK.R) and [results](test_around_branch_output.txt). However, it is worth noting that the four taxa subsets appear in a relative small fraction of the 852 UCE trees (114-183).
