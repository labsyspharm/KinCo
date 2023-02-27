library(shiny)
library(r3dmol)
library(bio3d)
library(colourpicker)
library(stringr)

kinco_meta <- read.table("KinCoMetaFile.csv",sep = ',',header = T)
kinases <- unique(kinco_meta$Symbol)


bio3d_docked<- function(homolog, poses){
  hom<- read.pdb(homolog)
  ps <- lapply(poses, function(x) read.pdb(x))
  docked <- cat.pdb(hom, ps[[1]], hom, rechain=FALSE, renumber=FALSE)
  if (length(ps)>1){
    
    
    for (i in 2:length(ps)) {
      
      docked <- cat.pdb(docked, ps[[i]], hom, rechain=FALSE, renumber=FALSE)
    }
  }
  return(docked)
}
