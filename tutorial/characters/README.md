TREE-QMC-BP Tutorial
---

1. Clone the repository.
```
git clone https://github.com/molloy-lab/TREE-QMC.git
```

2. Build wTREE-QMC and go to tutorial directory.
```
git clone https://github.com/molloy-lab/TREE-QMC.git
cd TREE-QMC/external/MQLib
make
cd ../..
g++ -std=c++11 -O2 \
    -I external/MQLib/include -I external/toms743 \
    -o treeqmc \
    src/*.cpp external/toms743/toms743.cpp \
    external/MQLib/bin/MQLib.a -lm \
    -DVERSION=\"$(cat version.txt)\"
```

3. Go to tutorial directory.
```
cd tutorial/characters
```

4. Run TREE-QMC-BP on the [example input data](4345ratites.nex). The nexus file contains presence/absence patterns for CR1 retrotransposon insertions from [Cloutier et al., *Syst Biol*, 2019](https://doi.org/10.1093/sysbio/syz019), plus the 44 additional characters found by [Simmons et al., *Mol Phy Evol*, 2022](https://doi.org/10.1016/j.ympev.2021.107344). The command for running TREE-QMC in `BP` mode is

```
../../treeqmc \
    --bp \
    --root galGal \
    --support \
    -i 4345ratites.nex \
    -o treeqmc-bp-4345ratites.tre
```
This command builds a species tree from the quartets induced by each character; see [Springer et al., *J Heredity*, 2020](https://doi.org/10.1093/jhered/esz076). If both `--support` and `--bp`  options are used, branch lengths will be computed using the MLE method described in [Molloy, Gatesy, Springer, *Syst Biol*, 2022](https://doi.org/10.1093/sysbio/syab086). 

**IMPORTANT: `BP` mode is for low-homoplasy bilallelic character matrices; however, you can run these analyses for multi-allelic characters (e.g., nucleotides) by replacing `--bp` by `--chars` in the command above.**

TREE-QMC-BP recovers the same species tree topology as [Cloutier et al., *Syst Biol*, 2019](https://doi.org/10.1093/sysbio/syz019), who estimated the species tree from gene trees using MP-EST (see the model tree below).

<p align="center">
<img src="model-species-tree-v2-annotated_low.jpg" alt="model" width="100%" height="auto"/>
</p>

The placement of Rhea in this species tree is debated. 
The quartet suport for this branch is `'q1=0.499331` (with `f1=18.650000` and `EN=37.350000`), which means that about half of characters with information about the placement of Rhea support making it sister to Kiwi+Emu+Cassowary. However, the low `EN` suggests that there is limited signal in the character matrix for resolving this branch.
The `--pcsonly` flag can be used to explore this futher by providing the quartet support for the placement of Rhea of **each** character.

5. The goal of **Partitioned Coalescence Support (PCS)**, described by [Gatesy et al., *Mol Phy Evol*, 2019](https://doi.org/10.1016/j.ympev.2019.106539), is to evaluate the quartet support for **each** character (or gene tree) in resolving a specific branch in the species tree. TREE-QMC can be used to compute PCS by providing an additional input: a species tree with `PCS` flagging the branch of interest; see [`species-tree-for-pcs.tre`](species-tree-for-pcs.tre) and [`unresolved-species-tree-for-pcs.tre`](unresolved-species-tree-for-pcs.tre). The command for running TREE-QMC in `PCS` mode is 

```
../../treeqmc \
    --bp \
    --pcsonly species-tree-for-pcs.tre \
    -i 4345ratites.nex \
    -o pcs-bp-4345ratites.tsv
```

The output shown is below (plus annotations (`#`) for each row with the tree it supports).
The four possible trees related to placement of Rhea are as follows:
* `t1 = A,B|C,D = Rhea,Kiwi+Emu+Cassowary|Tinamou,Chicken+Ostrich` (matches input species trees)
* `t2 = A,C|B,D = Rhea,Tinamou|Kiwi+Emu+Cassowary,Chicken+Ostrich`
* `t3 = Rhea,Chicken+Ostrich|Kiwi+Emu+Cassowary,Tinamou`
* `t4 = A,C,C,D = Rhea,Kiwi+Emu+Cassowary,Tinamou,Chicken+Ostrich` (no quartet / polytomy)

Thus, there are at most `|x| * |y| * |z| * |w| = 2 * 5 * 4 * 2 = 80` quartets with information about the placement of Rhea (although a character will induce fewer quartets if there are missing data).
Because each character contains can induce multiple quartets, it can split its vote between `t1`, `t2`, `t3`, or `t4`. Overall, the PCS analysis finds 61 characters (out of 4345 characters) that support at least one of `t1`, `t2`, and `t3`. If polytomies (`t4`) are ignored, 31 characters (51%) have quartets supporting only `t1`, 14 characters (23%) have quartets supporting only `t2`, and 16 characters (26%) have quartets supporting only `t3`. This is close to the estimates above for the branch (note `f1` above is computed by dividing column `f_xy|zw` by column `totalf` and then summing the result).

```
# x = rheAme,rhePen
# y = aptHaa,aptOwe,aptRow,casCas,droNov
# z = cryCin,tinGut,eudEle,notPer
# w = strCam,galGal
id	position	f_xy|zw	f_xz|yw	f_xw|yz	f_xyzw	totalf
3340	3341	30	0	0	0	30 # t1
3341	3342	40	0	0	0	40 # t1
3342	3343	40	0	0	0	40 # t1
3343	3344	40	0	0	0	40 # t1
3344	3345	30	0	0	0	30 # t1
4276	4277	24	0	0	16	40 # t1
4277	4278	24	0	0	16	40 # t1
4278	4279	24	0	0	16	40 # t1
4279	4280	18	0	0	12	30 # t1
4280	4281	12	0	0	8	20 # t1
4281	4282	16	0	0	24	40 # t1
4282	4283	8	0	0	12	20 # t1
4283	4284	16	0	0	24	40 # t1
4284	4285	12	0	0	18	30 # t1
4285	4286	16	0	0	24	40 # t1
4286	4287	16	0	0	24	40 # t1
4296	4297	0	4	0	6	10 # t2
4297	4298	0	0	24	16	40 # t3
4298	4299	0	0	16	24	40 # t3
4299	4300	0	40	0	0	40 # t2
4300	4301	0	10	0	10	20 # t2
4301	4302	0	0	30	0	30 # t3
4302	4303	0	0	30	0	30 # t3
4303	4304	0	0	40	0	40 # t3
4304	4305	0	0	30	0	30 # t3
4305	4306	0	0	40	0	40 # t3
4306	4307	0	0	40	0	40 # t3
4307	4308	40	0	0	0	40 # t1
4308	4309	40	0	0	0	40 # t1
4309	4310	0	40	0	0	40 # t2
4310	4311	0	40	0	0	40 # t2
4311	4312	0	0	16	24	40 # t3
4312	4313	0	0	8	12	20 # t3
4313	4314	0	0	12	18	30 # t3
4314	4315	0	0	16	24	40 # t3
4315	4316	0	0	12	18	30 # t3
4316	4317	0	0	24	16	40 # t3
4317	4318	0	0	18	12	30 # t3
4318	4319	0	0	18	12	30 # t3
4319	4320	24	0	0	16	40 # t1
4320	4321	24	0	0	16	40 # t1
4321	4322	24	0	0	16	40 # t1
4322	4323	24	0	0	16	40 # t1
4323	4324	24	0	0	16	40 # t1
4324	4325	18	0	0	12	30 # t1
4325	4326	16	0	0	24	40 # t1
4326	4327	16	0	0	24	40 # t1
4327	4328	12	0	0	18	30 # t1
4328	4329	16	0	0	24	40 # t1
4329	4330	8	0	0	24	32 # t1
4330	4331	16	0	0	24	40 # t1
4331	4332	16	0	0	24	40 # t1
4332	4333	0	24	0	16	40 # t2
4333	4334	0	18	0	12	30 # t2
4334	4335	0	16	0	24	40 # t2
4335	4336	0	12	0	18	30 # t2
4336	4337	0	16	0	24	40 # t2
4337	4338	0	12	0	18	30 # t2
4338	4339	0	16	0	24	40 # t2
4339	4340	0	12	0	18	30 # t2
4340	4341	0	16	0	24	40 # t2
```


