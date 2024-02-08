library(data.table)
#library(readxl)
library(tidyverse)
library(magrittr)
library(dplyr)

get_sample_info <- function(){
  return(fread("~/ipscAstrocyteNeuron/files/sample_map.tsv"))
}

get_sample_metadata <- function(){
  my_data <- read_excel("/home/kdsa/ipscAstrocyteNeuron/files/BX0175_extraction_data.xlsx")
  head(my_data,2)
  data <- my_data %>% 
    mutate(sample = gsub("_", "-", `Sample name`)) %>% dplyr::select(`sample` ,`Additional`, `Position`,`eRIN`) 
  data <- data.table(dplyr::rename(data, sample_desc = Additional))
  head(data, 2)
  ctrl <- get_ctrl_info()
  ctrl <- ctrl[, .(ctrl_id, ctrl_desc, Age_of_donor, Sex_of_donor)] 
  map <- fread("~/ipscAstrocyteNeuron/files/sample_map.tsv")
  map <- unique(map[,.(sample, indi)])
  map_ctrl <- merge(map, ctrl, by.x="indi", by.y="ctrl_id", all.x=T)
  rm(map, ctrl)
  head(map_ctrl,2)
  data <- merge(data, map_ctrl, by = "sample", all.x=T)
  
  return(data)
}

get_ctrl_info <- function(){
  my_data <- fread("/home/kdsa/ipscAstrocyteNeuron/files/from_MC/RNAseqcell layout_March2019_Minee_withID.tsv" )
  #head(my_data,2)
  return(my_data)
  
}
  

map_file2sample_info <- function(fname, sampleID_colname= "SampleID") {
  #the function will read the file , add in the sample info and return a data table contiangn the merged data
  #fname - file to which the sample info has to be added
  #sampleID_colname - the name of the column that contains the sampleID e.g. NM6286_E5-BX0175-037_S1
  
  #fname <- "~/ipscAstrocyteNeuron/results/from_multiqc/star_alignment.tsv"
  #sampleID_colname<-"Category"
  data <- fread(fname)
  setnames(data, old = sampleID_colname, new= "SampleID")
  
  data$id <- substr(data$SampleID, 1, 20)
  setnames(data, old = "SampleID" , new= sampleID_colname )
  
  head(data,2)
  
  #get the sample mapping
  sample_map <- get_sample_info()
  head(sample_map,2)
  data <- merge( data, sample_map, by.x= "id", by.y="sampleId_long")
  rm(sample_map)
  return(data)
}


map_sampleIds2sample_info <- function(sampleIds) {
  #the function will return sample info for the sampleIDs (e.g. NM6286_E5-BX0175-037_S1)
  #sampleIds <- colnames(txi$counts)
  data <- data.table(sampleIds_queried = sampleIds)
  #head(data,2)
  
  data$id <- substr(data$sampleIds_queried, 1, 20)
  
  head(data,2)
  
  #get the sample mapping
  sample_map <- get_sample_info()
  head(sample_map,2)
  data <- merge( data, sample_map, by.x= "id", by.y="sampleId_long")
  rm(sample_map)
  return(data)
}

#make_compact_names <- function(longVersion){
make_compact_names <- function(dt, colname){
  #longVersion <- comps
  l_dt <- copy(dt)
  setnames(l_dt, old=colname, new="t_long_version")
  l_dt$compact_version <-  gsub( "_1uM_syn", "", l_dt$t_long_version)
  l_dt$compact_version <-  gsub( "_treatment", "_treat", l_dt$compact_version)
  l_dt$compact_version <-  gsub( "cyte", "", l_dt$compact_version)
  setnames(l_dt, old = c("t_long_version","compact_version"), 
           new= c(colname, paste0(colname, "_compact")))
  return(l_dt)
}
