
#Script the examines the enrichment in the gene common to DEd and DGE
#run enrichment with the overlapping dge and Ded 
#1st – all genes tested for DGE (ii) all sig dge.  Single plot for astrocytes and coculture. 

#Based on teh decision to go with the bg of all dig DGe
#running the same with the DS - all DS bg and  query geneset as - sigDS adn editing


library(data.table)
library(here)
library(org.Hs.eg.db)
library(ReactomePA)
library(clusterProfiler)


library(UpSetR)
source(here("./R/common_dge.R"))
source(here("./R/common_clustProfiler_upsetR.R"))

in_path <- here("../results/mRNA_editing/dge_ds_overlap/")
out_res_clustProf <- paste0(in_path, "cluster_profiler/check4paper/" )

outPath_dge_go <- paste0(out_res_clustProf, "dge/go/")
outPath_dgeFC_go <- paste0(out_res_clustProf, "dgeFC/go/")


if(!dir.exists(outPath_dge_go)){
  dir.create(outPath_dge_go, recursive = T)
  #dir.create(outPath_dgeFC_go, recursive = T)
  #dir.create(outPath_reac)
}

outPath_ds_go <- paste0(out_res_clustProf, "ds/go/")
if(!dir.exists(outPath_ds_go)){
  dir.create(outPath_ds_go, recursive = T)
}


run_compareCluster_oraGo4ont <- function(gene_list2test, bg_genes, outfile_prefix, ont2check= c("BP", "CC"), plot_title = "" ){
  #compare cluster run for the genelist
  
  #compare cluster GO 4 ontology
  ret_ont <- lapply(ont2check, function(ontology){
    o_file <- paste0(outfile_prefix, ontology, ".png")
    ret <- run_compareCluster_go4ontology_genelist(gene_list2test, bg_genes, ontology,  o_file, title_suffix=plot_title)
  })
  return(T)
}

comparison_astr <- "astrocyte_treated_vs_untreated"
comparison_cocul <- "coculture_treated_vs_untreated"

#DGE 
dge_ed <-  fread(paste0(in_path, "dge_edit_comb_geneset.tsv"))
head(dge_ed, 2)
table(dge_ed$studied_in)
dge_studied <- dge_ed[studied_in %in% c("both", "dge"), ]
rm(dge_ed)
table(dge_studied$comparison)

bg_astr_genes <- dge_studied[ comparison == comparison_astr, gene_name]
bg_cocult_genes <- dge_studied[ comparison == comparison_cocul, gene_name]

reg_vals <- c("up", "down")

#### Sig DGE - no FC cut off applied
#######################################
#ora GO -ont

#comparing enrichment of genes with sig_in = both - in the astr vs coculture

#Background - taking all genes studied in  the astrocytes and co-culture
bg_bothcult_geneids <- union(dge_studied[ comparison == comparison_astr, geneid], dge_studied[ comparison == comparison_cocul, geneid]) 
bg_bothcult_genenames <- union(dge_studied[ comparison == comparison_astr, gene_name], dge_studied[ comparison == comparison_cocul, gene_name]) 


table( dge_studied$comparison)
q_astr_genes <- dge_studied[  comparison == comparison_astr & sig_in == "both" , gene_name] 
q_cocult_genes <- dge_studied[  comparison == comparison_cocul & sig_in == "both", gene_name] 

gene_list2test <- list("Astr-sig in DEd DGE"= q_astr_genes, 
                       "Co-culture-sig in DEd DGE" = q_cocult_genes)
bg_genes <- bg_bothcult_genenames #24470
outfile_prefix <- paste0( outPath_dge_go, "/bothCult_clust_cmp_sigin_bothDgeDed_enrichGO_")
run_compareCluster_oraGo4ont(gene_list2test, bg_genes, outfile_prefix, ont2check= c("BP", "CC") )
rm(gene_list2test, bg_genes, outfile_prefix)

#applying FC cut off of >= 2
q_astr_genes <- dge_studied[  comparison == comparison_astr & sig_in == "both" & abs(lfcShrink_log2FoldChange) >= 1 , gene_name] 
q_cocult_genes <- dge_studied[  comparison == comparison_cocul & sig_in == "both" & abs(lfcShrink_log2FoldChange) >= 1, gene_name] 

gene_list2test <- list("Astr-sig in DEd DGE"= q_astr_genes, 
                       "Co-culture-sig in DEd DGE" = q_cocult_genes)
bg_genes <- bg_bothcult_genenames

outfile_prefix <- paste0( outPath_dge_go, "/bothCult_clust_cmp_sigin_bothDgeDed_FCcutoff_enrichGO_")
run_compareCluster_oraGo4ont(gene_list2test, bg_genes, outfile_prefix, ont2check= c("BP", "CC") )
rm(gene_list2test, bg_genes, outfile_prefix)


#Background - taking sig DGE genes in both the astrocytes and co-culture
bg_bothcult_geneids <- union(dge_studied[ comparison == comparison_astr & reg %in% reg_vals, geneid], dge_studied[ comparison == comparison_cocul & reg %in% reg_vals, geneid]) 
bg_bothcult_genenames <- union(dge_studied[ comparison == comparison_astr & reg %in% reg_vals, gene_name], dge_studied[ comparison == comparison_cocul & reg %in% reg_vals, gene_name]) 

sig_dge_ed <- dge_studied[ geneid %in% bg_bothcult_geneids & sig_in =="both", ] 
table( sig_dge_ed$comparison)
length(unique(sig_dge_ed$gene_name)) 

q_astr_genes <- sig_dge_ed[  comparison == comparison_astr & sig_in == "both" , gene_name] 
q_cocult_genes <- sig_dge_ed[  comparison == comparison_cocul & sig_in == "both", gene_name] 

gene_list2test <- list("Astr-sig in DEd DGE"= q_astr_genes, 
                       "Co-culture-sig in DEd DGE" = q_cocult_genes)
bg_genes <- bg_bothcult_genenames
outfile_prefix <- paste0( outPath_dge_go, "/bothCult_clust_cmp_bgSigDGE_sigin_bothDgeDed_enrichGO_")
run_compareCluster_oraGo4ont(gene_list2test, bg_genes, outfile_prefix, ont2check= c("BP", "CC") )

rm(gene_list2test, bg_genes, outfile_prefix)


#applying FC cut off of >= 2
q_astr_genes <- sig_dge_ed[  comparison == comparison_astr & sig_in == "both" & abs(lfcShrink_log2FoldChange) >= 1 , gene_name] 
q_cocult_genes <- sig_dge_ed[  comparison == comparison_cocul & sig_in == "both" & abs(lfcShrink_log2FoldChange) >= 1, gene_name]

gene_list2test <- list("Astr-sig in DEd DGE"= q_astr_genes, 
                       "Co-culture-sig in DEd DGE" = q_cocult_genes)
bg_genes <- bg_bothcult_genenames

outfile_prefix <- paste0( outPath_dge_go, "/bothCult_clust_cmp_bgSigDGE_sigin_bothDgeDed_FCcutoff_enrichGO_")
run_compareCluster_oraGo4ont(gene_list2test, bg_genes, outfile_prefix, ont2check= c("BP", "CC") )
rm(gene_list2test, bg_genes, outfile_prefix)


#looking at the sig set in both analyses DGE DEd - enrichment of the genes in this catgory that are common to astr and cocult vs only in each culture
######################################################################################################################################################
#DGE & DEd studied in both
dge_ed <-  fread(paste0(in_path, "dge_edit_comb_geneset.tsv"))
head(dge_ed, 2)

dge_ed_studied_both <- dge_ed[studied_in == "both", ]
rm(dge_ed)

#No FC applied -   done previously and rsults in prev run - /sigInBothDgeDed_clust_cmp_genesCommon2astrCocul_onlyAstr_onlyCocul_*

#Applying FC
#comparing the genes intersection of the sig in both in the astr and the cocul vs those that are not
sigInBoth_astr <- dge_ed_studied_both[  comparison == comparison_astr & sig_in == "both" & abs(lfcShrink_log2FoldChange) >= 1, gene_name] 
sigInBoth_cocul <-  dge_ed_studied_both[  comparison == comparison_cocul & sig_in == "both" & abs(lfcShrink_log2FoldChange) >= 1, gene_name] 
sigInBoth_comm2astrCocul <- intersect(sigInBoth_astr, sigInBoth_cocul)
sigInBoth_onlyInAstr <- sigInBoth_astr[ ! sigInBoth_astr %in% sigInBoth_comm2astrCocul ]
sigInBoth_onlyInCocul <- sigInBoth_cocul[ ! sigInBoth_cocul %in% sigInBoth_comm2astrCocul ]

#comparing the genes sig in both dge and Ded - common to the astr and coculture, only in astr and only in cocult
gene_list2test <- list("Sig DGE DEd - in Astro & Cocul " = sigInBoth_comm2astrCocul, 
                       "Sig DGE DEd - in Astro only" = sigInBoth_onlyInAstr, 
                       "Sig DGE DEd - in Cocul only" = sigInBoth_onlyInCocul)
bg_genes <- unique(dge_ed_studied_both$gene_name) #16414
outfile_prefix <- paste0( outPath_dge_go, "/sigInBothDgeDed_clust_cmp_genesCommon2astrCocul_onlyAstr_onlyCocul_FCcutoff_enrichGO_")
run_compareCluster_oraGo4ont(gene_list2test, bg_genes, outfile_prefix, ont2check= c("BP", "CC") )


#chnaging bg genes to those significant ( union of astr and cocult)
bg_genes <- unique(dge_ed_studied_both[ reg %in% reg_vals, gene_name])

outfile_prefix <- paste0( outPath_dge_go, "/sigInBothDgeDed_clust_cmp_bgSigDGE_genesCommon2astrCocul_onlyAstr_onlyCocul_FCcutoff_enrichGO_")
run_compareCluster_oraGo4ont(gene_list2test, bg_genes, outfile_prefix, ont2check= c("BP", "CC") )


#DS
######

#of the genes with DS what is special about the ones edited
ds_ed <-  fread(paste0(in_path, "ds_edit_comb_geneset.tsv"))
ds_ed_comm <- ds_ed[studied_in == "both", ]

#Background - taking sig DS genes in both the astrocytes and co-culture
bg_bothcult_genenames <- union(ds_ed_comm[ comparison == comparison_astr & ds_sig =="T", genes], ds_ed_comm[ comparison == comparison_cocul & ds_sig =="T", genes]) 
length(unique(ds_ed_comm[ ds_sig == "T", genes])) 


q_astr_genes <- ds_ed_comm[  comparison == comparison_astr & sig_in == "both" , genes] #95
q_cocult_genes <- ds_ed_comm[  comparison == comparison_cocul & sig_in == "both", genes] #81

gene_list2test <- list("Astr-sig in DEd DS"= q_astr_genes, 
                       "Co-culture-sig in DEd DS" = q_cocult_genes)
bg_genes <- bg_bothcult_genenames
outfile_prefix <- paste0( outPath_ds_go, "/bothCult_clust_cmp_bgSigDS_sigin_bothDSDed_enrichGO_")
run_compareCluster_oraGo4ont(gene_list2test, bg_genes, outfile_prefix, ont2check= c("BP") )
# No enrichment found in any of gene cluster, please check your input... 
run_compareCluster_oraGo4ont(gene_list2test, bg_genes, outfile_prefix, ont2check= c("CC") )
# No enrichment found in any of gene cluster, please check your input...
rm(gene_list2test, bg_genes, outfile_prefix)

