Quartet Input Tutorial
======================

This tutorial shows how to run TREE-QMC to reconstruct a tree from input quartets.

---

If you have not already, build TREE-QMC. Also, be sure that the TREE-QMC excutable is on your path; see [tips](../README.md). 

To begin the tutorial, go to directory containing the example input data by typing
```
cd TREE-QMC/tutorial/quartets
```
TREE-QMC can take quartets in QMC, QFM, or PhyloNetworks qCF format as input.


QMC format
---
 To see the QMC format, type
```
head -n3 qmc.txt 
```
This command should return
```
leaf1,leaf10|leaf2,leaf3:10.000000
leaf1,leaf2|leaf10,leaf3:78.000000
leaf1,leaf3|leaf10,leaf2:12.000000
```
Now type
```
tree-qmc -i qmc.txt --quartets --quartetformat "___,___|___,___:___"
```
to build a tree. 

**Note:** The above command uses the `n2` artificial taxon normalization scheme by default. To change it to `n0`, add the flag `--norm_atax 0`. To learn more about normalization, see [Han & Molloy, *Genome Res*, 2023](https://doi.org/10.1101/gr.277629.122).


QFM format
---
To see the QFM format, type
```
head -n3 qfm.txt 
```
This command should return
```
((leaf1,leaf10),(leaf2,leaf3));10.000000
((leaf1,leaf2),(leaf10,leaf3));78.000000
((leaf1,leaf3),(leaf10,leaf2));12.000000
```
Now type 
```
tree-qmc -i qfm.txt --quartets
```
to build a tree.

PhyloNetworks qCF format
---
To see the PhyloNetworks qCF format, type
```
head -n2 tableCF.txt
```
This command should return 
```
t1,t2,t3,t4,CF12_34,CF13_24,CF14_23,ngenes
leaf1,leaf10,leaf2,leaf3,0.1,0.78,0.12,100.0
```
Now type 
```
tree-qmc -i tableCF.txt --quartets
```
to build a tree.
