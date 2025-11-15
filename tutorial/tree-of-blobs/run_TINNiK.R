# https://www.rdocumentation.org/packages/MSCquartets/versions/3.2

library('MSCquartets')

# Run TINNiK
args <- commandArgs(trailingOnly=TRUE)
gene_trees <- read.tree(file=args[1])
taxon_names <- taxonNames(gene_trees)
output <- TINNIK(gene_trees,
                 alpha=as.numeric(args[3]),
                 beta=as.numeric(args[4]),
                 plot=FALSE)
write.tree(output$ToB, file=args[2])
