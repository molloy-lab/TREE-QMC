#!/bin/bash
declare -a gene_tree=(
  "avian_uce_trees_3679"
  "1_best_cds5756"
  "butterfly_bs"
  "zhang_cetacean"
  "1_best_intron4871"
  "avian14k_bs"
  "CNEE_abayes_gene_trees_sorted"
  "zhang_papilionidae"
  "1kp-AA-genes-bayes"
  "swallowtail_bayes"
  "1kp-AA-genes-bs"
  "bee_bs"
  "whale_bayes"
)

g++ -std=c++17 -O3 -I external/MQLib/include -I external/toms743 -I external/parlaylib/include -o tree-qmc -pthread -mcx16 -march=native src/*.cpp external/toms743/toms743.cpp external/MQLib/bin/MQLib.a -lm -DVERSION=\"$(cat version.txt)\"

for input in "${gene_tree[@]}"; do
  echo Running on ${input}.tre
  ./tree-qmc --hybrid --bootstrap -i /ssd1/treeqmc/${input}.tre -o ${input}.out
  echo
done