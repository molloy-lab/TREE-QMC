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

## ----eval=FALSE---------------------------------------------------------------
# # read text file of gene trees and count quartets on them
# 
# gts<-read.tree(file = 'genetreefile')
# tableLeopardusLescroart=quartetTable(gts)

## -----------------------------------------------------------------------------
# load data file containing quartet counts for Leopardus data set supplied with MSCquartets package
# These counts are will be accessed as `tableLeopardusLescroart`.

data(tableLeopardusLescroart)

## -----------------------------------------------------------------------------
# perform TINNIK analysis for gene trees, using defaults

output<-TINNIK(tableLeopardusLescroart)

# save table of quartet information with p-values

pT<-output$pTable

## -----------------------------------------------------------------------------
# perform improved TINNIK analysis to infer the tree of blobs

output<-TINNIK(pT, alpha=5e-29,beta = 0.95)

## -----------------------------------------------------------------------------
# run TINNIK to infer the tree of blobs

output<-TINNIK(pT, alpha=5e-29,beta = 0.95)

## -----------------------------------------------------------------------------
# rename output 

pT<-output$pTable #quartet count data with p-values for tests
ToB<-output$ToB   #the TINNIK tree of blobs

## -----------------------------------------------------------------------------
# perform NANUQ analysis for table of quartet information and p-values

D<-NANUQ(pT, alpha = 5e-29,beta = 0.95) # Run the NANUQ routine
NN<-neighborNet(D$dist) # Run the NeighborNet algorithm on the NANUQ distance
plot(NN) # plot the splits-graph with neighborNet

## -----------------------------------------------------------------------------
  #Label internal nodes of the tree of blobs, and plot

  ToB<-labelIntNodes(ToB)

## -----------------------------------------------------------------------------
  # resolve node 18

  resC18<-resolveCycle(ToB,18,pT,alpha=5e-29,beta=0.95)

## -----------------------------------------------------------------------------
  # resolve node 2O

  resC20<-resolveCycle(ToB,20,pT,alpha=5e-29,beta=0.95)

## -----------------------------------------------------------------------------
#Fully resolve the tree of blobs to a level-1 network

resN<-resolveLevel1(ToB=output$ToB, pTable=output$pTable, alpha=5e-29, beta=0.95, 
                    distance="NANUQ")

