#generating the data for the ADAR - alternative splciing ggtranscript plot
#code from create_se_obj_4sashimi_plots.R
library(here)
source(here("./R/common_sample.R"))
library(data.table)

outPath <- here("../results/paper/data4plots/")
gene_of_interest <- "ADAR"  #goi =gene of interest

leafcutter_res_path <- here("../results/leafcutter_additionalFilters_postdeconv/")
dt_cluster_ann <- fread(paste0(leafcutter_res_path, "/diff_splicing/combGrp_clust_junct_anno_singleGene.tsv")) #cluster info and anotation

head(dt_cluster_ann, 2)
dt_cluster_ann <- dt_cluster_ann[ genes == gene_of_interest & p.adjust < 0.05, ]

dt_int <- fread(paste0(leafcutter_res_path, "/diff_splicing/combGrp_clust_effSize_singleGene.tsv"))
dt_int <- dt_int[ p.adjust < 0.05, ] # as it has the cluster info - keeping intron info in signinficant clusters

#dt_int <- dt_int[ p.adjust < 0.05 &  comparison == comp , ] # as it has the cluster info - keeping intron info in signinficant clusters
dt_int <- dt_int[  intron %in% dt_cluster_ann$intron, ] #keeping introns where junction has required annotation
head(dt_int, 2)
#add the junction_id
setDT(dt_int)[dt_cluster_ann, junction_id := junction_id, on = c(intron="intron") ]
dt_int <- dt_int[, c("strand") := tstrsplit(intron, "_", fixed=TRUE, keep = c(3))]
dt_int[ ,chr := gsub("chr", "", chr)]
setnames(dt_int, old = "1uM_syn_oligomer", new = "treated")
setnames(dt_int, old = "no_treatment", new = "untreated")
cols2Sel <- c( "junction_id", "chr", "start", "end", "strand", "treated", "untreated", "comparison")

intr_usage <- dt_int[ , .SD, .SDcols= cols2Sel ]
rm(dt_cluster_ann, dt_int)
#
id_cols <- c( "junction_id", "chr", "start", "end", "strand", "comparison")
melted_intr_usage <- melt(intr_usage, id=id_cols )
setnames(melted_intr_usage, old = c("variable", "value"), new = c("treatment", "intron_usage"))
melted_intr_usage$culture <- gsub("_treatment_vs_no_treatment", "", melted_intr_usage$comparison)
melted_intr_usage[ culture == "astrocyte_neuron", culture := "coculture"]
melted_intr_usage$treatment <- paste0(melted_intr_usage$culture, "_", melted_intr_usage$treatment)

#get the transcripts of overlap with the junctions
junct_overlap_path <- paste0(leafcutter_res_path , "/diff_splicing/junct_of_interest")
junct_overlap_files <- list.files(junct_overlap_path, pattern = "txOverlapJunc.tsv", full.names = T)
junct_overlap_files <- junct_overlap_files[ grep(gene_of_interest, junct_overlap_files)]
dt_juncOverlapTx <- lapply(junct_overlap_files, function(f){
  #f <- junct_overlap_files[1]
  l_dt <- fread(f)
  l_dt$comparison <- gsub("_treat_vs_no_treat_txOverlapJunc.tsv" , "", tstrsplit( basename(f), split = "-", keep=2) )
  return(l_dt)
})
dt_juncOverlapTx <- rbindlist(dt_juncOverlapTx)
dt_juncOverlapTx <- dt_juncOverlapTx[ transcript_biotype == "protein_coding", ]
dt_juncOverlapTx[ comparison == "_astr_neuron", comparison:="astrocyte_neuron_treatment_vs_no_treatment"]
dt_juncOverlapTx[ comparison == "_astr", comparison:="astrocyte_treatment_vs_no_treatment"]
dt_juncOverlapTx <- dt_juncOverlapTx[ , .(external_transcript_name, junction_id, comparison)]
melted_intr_usage <- merge(melted_intr_usage, dt_juncOverlapTx, by = c("comparison", "junction_id")) #significant junctions
setnames(melted_intr_usage, old = "external_transcript_name", new= "transcript_name")
rm(dt_juncOverlapTx, junct_overlap_files)
#save the jucntion info
fwrite(melted_intr_usage, paste0(outPath, "adar_junctionTx_data.tsv"), quote = F, sep="\t")