library(reticulate)
library(dplyr)
library(stringr)
#use_python("/opt/anaconda3/bin/python3.7") #You may or may not need to change the python source based on your computer read https://rstudio.github.io/reticulate/ 
#pd <- import('pandas')
#Change the compound-kinase pair
#tgtgeneid = '27'
#lgdinchikey = 'ZBNZXTGUTAYRHI-UHFFFAOYSA-N'
#modelname = 'ABL1_HUMAN_D0_1OPL_A'
#poseID = '0-18'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


unzip_lgd <- function(lgdtmp){
  if (!dir.exists(lgdtmp)){
    dir.create(lgdtmp)
  }
  print(lgdtmp)
  unzip(zipfile = paste0(lgdtmp,'.zip'), exdir = lgdtmp)
  return(paste0(lgdtmp,'/'))
}

export_PDBQT <- function(refpdbqt, coords, new_file){
  lines = readLines(refpdbqt)
  acc = 1
  
  sink(new_file)
  for (aline in lines){
    if (grepl('ATOM', x = aline) || grepl('UNL', x= aline)){
      new_line <- aline
      new_coord = unlist(lapply(as.numeric(coords[acc,]), as.character))
      new_coord = str_pad(new_coord,width=8, side = 'left')
      substr(new_line, 30, 37) <- new_coord[1]
      substr(new_line, 38, 45) <- new_coord[2]
      substr(new_line, 46, 53) <- new_coord[3]
      acc = acc +1
    } else{
      new_line <- aline
    }
    cat(new_line,sep="\n")
  }
  
  sink()
}

#posedf <- data.frame()

# datasetdir = './' 
# pairpath = paste0(datasetdir,tgtgeneid,'/',lgdinchikey)
# output_dir = paste0(datasetdir,'/exported_structures/')
# if (!dir.exists(output_dir)){
#   dir.create(output_dir)
# }
# 
# pairdir <- unzip_lgd(pairpath)
# pdbqt = paste0(pairdir,'/PDBQT/',lgdinchikey,'.pdbqt')
# posedf = pd$read_pickle(paste0(pairdir,'pose.df'))
# pdbqt = paste0(pairdir,'/PDBQT/',lgdinchikey,'.pdbqt') #read pose coordinate
# if (!dir.exists(paste0(output_dir,tgtgeneid,'-',lgdinchikey))){
#   dir.create(paste0(output_dir, tgtgeneid, '-', lgdinchikey))
# }
# 
# ligandchain = posedf %>% filter(ModelName == modelname, PoseID ==poseID) %>%.$LigandChain%>%.[[1]]
# ligandchain = data.frame(x = ligandchain$x, y = ligandchain$y, z = ligandchain$z)
# 
# posepath = paste0(output_dir, tgtgeneid, '-', lgdinchikey, '/',modelname,'-',poseID,'_lgd.pdb')
#export_PDBQT(pdbqt, ligandchain, posepath)



