source("~/ipscAstrocyteNeuron/ipscAstroNeu/R/common_sample.R")

getTxInfo_ensembl <- function(enstIDs, versionNo) {
  #filter applioed to transcript_id_version eg ENST00000571311.5
  
  #returns the info using ensembl 93
  #enstIDs <- head(rownames( txi$counts))
  
  library(biomaRt)
  
  options(ucscChromosomeNames =FALSE)
  # For Ensembl 93: (GRCh38) -
  
  ens <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", version =versionNo)
  #for a gene name
  geneDetails <- getBM(attributes=c("ensembl_gene_id","external_gene_name",
                                    "ensembl_transcript_id", "ensembl_transcript_id_version", 
                                    "gene_biotype", "transcript_biotype", 
                                    "transcript_start", "transcript_end", "chromosome_name"),
                       filters="ensembl_transcript_id_version",values = enstIDs, mart=ens, uniqueRows = T)
  
  
  return(geneDetails)
}

getTxInfo_ensembl_4txid <- function(enstIDs, versionNo=93) {
  
  #Copy of getTxInfo_ensembl, to get additonal tx info AND
  ##Filter applied to transcript_id e.g. ENST00000571311 
  
  #returns the gene info using ensembl 93 default 
  # and transcript name , tsl if  txname is T
  #enstIDs <- intron_in_tx # junctions
  library(biomaRt)
  
  options(ucscChromosomeNames =FALSE)
  # For Ensembl 93: (GRCh38) -
  
  ens <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", version =versionNo)
  geneDetails <- getBM(attributes=c("ensembl_gene_id","external_gene_name",
                                    "gene_biotype",
                                    "ensembl_transcript_id", "ensembl_transcript_id_version", 
                                     "transcript_biotype", 
                                    "external_transcript_name", "transcript_tsl",
                                    "transcript_start", "transcript_end", "chromosome_name", "ensembl_peptide_id" ),
                       filters="ensembl_transcript_id",values = enstIDs, mart=ens, uniqueRows = T)
    
  return(geneDetails)
}



#get the quant file name with path - function from the nucCyt project
#returns a data frame with the file name( incl path) and sample
getSalmonQuantFiles <- function(path2quant){
  
  sFiles <- data.frame(fname= list.files(path2quant, pattern="quant.sf", recursive = T, full.names = T),
                       sample=basename(dirname(list.files(path2quant, pattern="quant.sf", recursive = T, full.names = T))), stringsAsFactors = F)
  
  return(sFiles)
}

getSampleInfo4DeseqColdata <- function(sampleIds){
   sampleInfo <-map_sampleIds2sample_info(sampleIds)
   head(sampleInfo,1)
   data <- sampleInfo[, .(sampleIds_queried, cells, treatment, indi,  replicate, sample)]
   mData <- data.table(get_sample_metadata())
   head(data,2)
   head(mData,2)
   #setDT(data)[mData,Position := Position, on=(sample = "sample") ]
   setDT(data)[mData,eRIN := eRIN, on=(sample = "sample") ]
   setDT(data)[mData,age_of_donor := Age_of_donor, on=(sample = "sample") ]
   setDT(data)[mData,sex_of_donor := Sex_of_donor, on=(sample = "sample") ]
   
   data$sample <- NULL
   #sQ <- data$sampleIds_queried
   #stopifnot(identical(sQ, sampleIds) ) # want the sampleids in the same order
   
   setDF(data)
   rownames(data) <- data$sampleIds_queried
   data$sampleIds_queried <- NULL
   
   rm(sampleInfo)
   return(data)
   

}
