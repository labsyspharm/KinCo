library(bio3d)
library(r3dmol)
library(reticulate)
library(testthat)
#library(iterators)
library(rlist)
source("molViewer.R")

#use_python("/opt/anaconda3/bin/python3.7", required = T)
stopifnot(py_module_available("Bio")==TRUE) #should be true

Bio <- import("Bio")
#BioPDB <- import('Bio.PDB')
PDB <- Bio$PDB
#SeqIO <- import('Bio.SeqIO')
#SeqIO <- Bio$SeqIO
#PDBIO <- import("Bio.PDB.PDBIO.PDBIO")
PDBIO <- Bio$PDB$PDBIO
PDB_parser <- PDB$PDBParser()
homolog_1 <- PDB_parser$get_structure("homolog_1", "data/PAK1_HUMAN_D0_4DAW_A_tgt.pdb")
homolog_2 <- PDB_parser$get_structure("homolog_2", "data/STK3_HUMAN_D0_4LG4_F_tgt.pdb")

#io <- PDBIO()
#io$set_structure(homolog_1)
#io$save("homolog_1.pdb")


#get atoms list
atoms1 <- homolog_1$get_atoms()
atoms2 <- homolog_2$get_atoms()

latoms1 <- c()
latoms2 <-c()


for (a in iterate(atoms1)) {
  #list.append(latoms1, a)
  latoms1 <- c(latoms1, a)
}
  
for (a in iterate(atoms2)){
  latoms2 <- c(latoms2, a)
}
  

super_imposer <- PDB$Superimposer()
super_imposer$set_atoms(fixed = latoms1, moving = latoms2)
super_imposer$apply(latoms2)

io <- PDBIO()
io$set_structure(homolog_2)
io$save("homolog_2.pdb")

io <- PDBIO()
io$set_structure(homolog_1)
io$save("homolog_1.pdb")
