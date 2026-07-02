TOB-QMC Tutorial
=================

This tutorial shows how to run NetCS to reconstruct a level-1 semi-directed phylogenetic network by (partially)resolving blobs of degree at least 5. NetCS is built within TOB-QMC.

BUILD
=====

Requirements:
* [cmake](https://cmake.org/download/)
* [R](https://www.r-project.org/)
* [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)
* [RInside](https://cran.r-project.org/web/packages/RInside/index.html) 
* [MSCquartets](https://cran.r-project.org/web/packages/MSCquartets/)

**Step 1.** If cmake version 3.18 or later isn't already installed on your system, install cmake by following the instructions on the [cmake official website](https://cmake.org/download/). If you are using a cluster, you may be able to load cmake as a module. Type 
```
module avail cmake
```
to see the available modules. If cmake is available, load it by typing
```
module load <module name>
```

**Step 2.** If R isn't already installed on your system, install R following the instructions on the [R official website](https://www.r-project.org/). If you are using a cluster, you may be able to load R as a module. Type `module avail R` to see the available packages. If R is available, load it by typing `module load <module name>`.

**Step 3.** Start the R console by typing `R` into the terminal as a commandline instruction and install the depedencies with the following R commands:
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

**Step 4.** Now you are ready to download and build NetCS. Just type the following commands:
```
git clone https://github.com/molloy-lab/TREE-QMC.git
cd TREE-QMC
mkdir -p build && cd build
cmake ..
cmake --build . -j
```

**Step 5.** Lastly, add `build` directory to your `$PATH` environment variable by typing
```
export PATH="$(pwd):$PATH"
cd ../..
```
so that your system can find `tree-qmc` from other directories (and then return to the next directory up so you are ready for the tutorial below). Even better, add the following line
```
export PATH="<path to TREE-QMC>:$PATH"
```
to your `~/.bash_profile` file or your `~/.bashrc` file so that your system can find `tree-qmc` whenever you start a new terminal instance.


TUTORIAL
========
In this tutorial, we will resolve all blobs with degree of at least 5 in a tree-of-blobs.

To begin, go to tutorial directory.
```
cd TREE-QMC/tutorial/NetCS
```

NetCS has two input files: (1)Input gene trees, and (2) Input tree of blobs and output the reconsturcted level-1 semi-directed netowrk 

In this example, ```gt.trees``` contains the gene trees, ```tob.tree`` contains the input tree of blobs. 

To run NetCS, type: 
```
tree-qmc \
    --network \
    -i gt.trees \
    --at tob.tree \
    -o output_net.nwk
```