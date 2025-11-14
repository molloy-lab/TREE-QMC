# https://www.rdocumentation.org/packages/MSCquartets/versions/3.2

library('MSCquartets')

# Test for hybridization around branch of interest
gene_trees <- read.tree("nomiinae_gene_trees.tre")

rotate_taxon_set <- c("Lipotriches_collaris",
	                  "Lipotriches_justiciae",
	             	  "Afronomia_circumnitens",
	             	  "Macronomia_clavisetis",
	             	  "Dieunomia_heteropoda",
	             	  "Dieunomia_triangulifera",
	             	  "Pachynomia_tshibindica",
	             	  "Pachynomia_amoenula",
	             	  "Pseudapis_oxybeloides",
	             	  "Pseudapis_siamensis",
	             	  "Pseudapis_riftensis",
	            	  "Pseudapis_pandeana",
	            	  "Pseudapis_kenyensis",
	            	  "Pseudapis_cinerea",
	            	  "Steganomus_junodi",
	             	  "Hoplonomia_elliotii",
	            	  "Austronomia_australica",
	            	  "Curvinomia_chalybeata",
	            	  "Acunomia_melanderi",
	           	      "Lasioglossum_albipes",
	             	  "Nomiapis_diversipes",
	             	  "Nomiapis_bispinosa")

for (taxon in rotate_taxon_set) {
	four_taxon_set <- c(taxon,
		                "Stictonomia_sangaensis",
		                "Stictonomia_aliceae",
		                "Stictonomia_schubotzi")

	qCFs <- quartetTable(gene_trees, four_taxon_set)
	test <- quartetTreeTest(qCFs[5:7], "T3")
	print(four_taxon_set)
	print(qCFs[5:7])
	print(qCFs[5] + qCFs[6] + qCFs[7])
	print(test$p.value)
	print("")
}

