## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>")
  
knitr::opts_chunk$set(fig.align = "center", 
                      fig.show = "hold",
                      out.width = "55%",
                      fig.width = 7,  
                      fig.height = 6)



## ----setup--------------------------------------------------------------------
library(MSCquartets)

## -----------------------------------------------------------------------------
# read text file of gene trees supplied with MSCquartets package

gts=read.tree(file = system.file("extdata","dataPapioniniVanderpool",package="MSCquartets"))


## -----------------------------------------------------------------------------

# perform initial TINNIK analysis for gene trees, using defaults

output=TINNIK(gts)

# save table of quartet information and p-values

pT=output$pTable


## -----------------------------------------------------------------------------

TINNIK(pT, alpha=.05, beta=1e-40)


## -----------------------------------------------------------------------------

TINNIK(pT, alpha=.01, beta=.95)


## -----------------------------------------------------------------------------

TINNIK(pT, alpha=.02, beta=.95)


