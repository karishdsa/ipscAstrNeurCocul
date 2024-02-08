
library(data.table)
library(here)
library(org.Hs.eg.db)
library(ReactomePA)
library(clusterProfiler)
library(UpSetR)
source(here("./R/common_dge.R"))
source(here("./R/common_clustProfiler_upsetR.R"))
out_res_clustProf <- here( "../results/cluster_profiler/" )
out_dge_path <- here("../results/dge_postdeconv/")

outPath_go <- paste0(out_res_clustProf, "dge/go/")
outPath_reac <- paste0(out_res_clustProf, "dge/reactome/")

if(!dir.exists(outPath_go)){
  dir.create(outPath_go, recursive = T)
  dir.create(outPath_reac)
 
}


#DGE

dge_all <- fread( paste0(out_dge_path, "correct_PC234/dge_astr_cocult_all.tsv"))
#head(dge_all, 2)
#dge_genes <- dge_all[ padj < 0.05, ]

comparison_astr <- "astrocyte_treated_vs_untreated"
comparison_cocul <- "coculture_treated_vs_untreated"

nrow(dge_all) #46883
dge_all <- dge_all[ , .SD[which.min(padj)], by = c("comparison", "gene_name")] #th
nrow(dge_all) #46814

#FC of >=2
fc2_abs_LFC_cutoff <- 1

#getting the genes at FDR <5% and FC>=2
dge_astr_all <- dge_all[ comparison== comparison_astr, ]
dge_cocult_all <- dge_all[ comparison== comparison_cocul, ]

#apply min padj


bg_astr_genes <- dge_astr_all$gene_name
bg_cocult_genes <- dge_cocult_all$gene_name

q_astr_genes <- dge_astr_all[ padj < 0.05 &
                              abs(lfcShrink_log2FoldChange) >= fc2_abs_LFC_cutoff, gene_name]
q_cocult_genes <- dge_cocult_all[ padj < 0.05 &
                               abs(lfcShrink_log2FoldChange) >= fc2_abs_LFC_cutoff, gene_name]

#GO enrichment analysis
#########################
##Gsea   #all genes tested
head(dge_astr_all, 2)
original_gene_list <- dge_astr_all[ ,lfcShrink_log2FoldChange]
names(original_gene_list) <- dge_astr_all$gene_name

run_gsea4genelist(original_gene_list,  paste0(outPath_go, "astr_gse_go.png") , n_show_cate = 30 ) 

original_gene_list <- dge_cocult_all[ ,lfcShrink_log2FoldChange]
names(original_gene_list) <- dge_cocult_all$gene_name
run_gsea4genelist(original_gene_list,  paste0(outPath_go, "cocul_gse_go.png") , n_show_cate = 30 )



#ORA
outfile_pref <- paste0( outPath_go,"/astr_ora")

#Running for all ont - previously giving an error. Therefore exlcuding it. 
run_go_ora_exclAll(q_astr_genes, bg_astr_genes, outfile_pref, n_show_category=30,  png_width =10, png_ht=10, png_units="in", plot_font_size =12, png_res=300) #runs the enrichGO() for each ontology folloewd by simplify()
# function also runs enrichGO() for 'ALL' ontology - simplify() does not work on this


outfile_pref <- paste0( outPath_go,"/cocul_ora")
run_go_ora_exclAll(q_cocult_genes, bg_cocult_genes, outfile_pref, n_show_category=30,  png_width =10, png_ht=10, png_units="in", plot_font_size =12, png_res=300) #runs the enrichGO() for each ontology folloewd by simplify()

#Running the enrich Reactome
############################
outfile_pref <- paste0( outPath_reac,"/astr_ora")
run_reac_ora (q_astr_genes, bg_astr_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300)

outfile_pref <- paste0( outPath_reac,"/cocul_ora")
run_reac_ora (q_cocult_genes, bg_cocult_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300)


#GSEA - reac
l_dge_astr_all <- copy(dge_astr_all)
setnames(l_dge_astr_all, old= c("gene_name", "lfcShrink_log2FoldChange" ), new=c("gene_symbol", "avg_log2FC"))  
names(l_dge_astr_all)
run_gsea_reac(l_dge_astr_all, paste0(outPath_reac, "astr_gse_reac.png") , n_show_cate = 30, png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300 )

rm(l_dge_astr_all)
l_dge_cocul_all <- copy(dge_cocult_all)
setnames(l_dge_cocul_all, old= c("gene_name", "lfcShrink_log2FoldChange" ), new=c("gene_symbol", "avg_log2FC"))  
names(l_dge_cocul_all)
run_gsea_reac(l_dge_cocul_all, paste0(outPath_reac, "cocul_gse_reac.png") , n_show_cate = 30, png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300 )
rm(l_dge_cocul_all)


#up and down regulated
#gsea


#ora & compCluster()
run_reac_ora_compclust_4up_dn_reg_genes <- function( dge_dt, abs_lfc_cutoff, out_path_prefix, cmpclust_listnameprefix){
  #runs ora and compareCluster() for reactome pathways
  #out_path_prefix <- paste0( outPath_reac,"/astr")
  #cmpclust_listnameprefix <- "Astrocytes"
  q_up <- dge_dt[ padj < 0.05 &
                               abs(lfcShrink_log2FoldChange) >= abs_lfc_cutoff &
                               reg == "up", ]
  q_dn <- dge_dt[ padj < 0.05 &
                               abs(lfcShrink_log2FoldChange) >= abs_lfc_cutoff &
                               reg == "down", ]
  
  bg_genes <- dge_dt$gene_name
  
  outfile_pref <- paste0( out_path_prefix,"_ora_up_reg")
  title_suffix <- "- Up regulated genes"
  run_reac_ora (q_up$gene_name, bg_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300, title_suffix)
  
  outfile_pref <- paste0( out_path_prefix,"_ora_dn_reg")
  title_suffix <- "- Down regulated genes"
  run_reac_ora (q_dn$gene_name, bg_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300, title_suffix)
  

  ofile_png <- paste0( out_path_prefix, "_clust_cmp_enrichPath_up_down_reg.png")
  run_compareCluster_reac(q_up$gene_name, qset1_listname= paste0(cmpclust_listnameprefix, "- Up regulated"), q_dn$gene_name, qset2_listname =paste0(cmpclust_listnameprefix, "- Down regulated"), bg_genes, ofile_png)

}


outfile_pref <- paste0( outPath_reac,"/astr")
run_reac_ora_compclust_4up_dn_reg_genes( dge_astr_all, fc2_abs_LFC_cutoff, outfile_pref, cmpclust_listnameprefix = "Astrocytes")

outfile_pref <- paste0( outPath_reac,"/cocul")
run_reac_ora_compclust_4up_dn_reg_genes( dge_cocult_all, fc2_abs_LFC_cutoff, outfile_pref, cmpclust_listnameprefix = "Coculture")

#running compareCluster() for the astr and cocult up and down togehter.

qastr_up <- dge_astr_all[ padj < 0.05 &
                  abs(lfcShrink_log2FoldChange) >= fc2_abs_LFC_cutoff &
                  reg == "up", .(gene_name)]
qastr_dn <- dge_astr_all[ padj < 0.05 &
                  abs(lfcShrink_log2FoldChange) >= fc2_abs_LFC_cutoff &
                  reg == "down", .(gene_name) ]

qcocul_up <- dge_cocult_all[ padj < 0.05 &
                  abs(lfcShrink_log2FoldChange) >= fc2_abs_LFC_cutoff &
                  reg == "up",  .(gene_name)]
qcocul_dn <- dge_cocult_all[ padj < 0.05 &
                  abs(lfcShrink_log2FoldChange) >= fc2_abs_LFC_cutoff &
                  reg == "down",  .(gene_name)]

bg_genes <- unique(dge_astr_all$gene_name, dge_cocult_all$gene_name) 

qastr_up <- bitr(qastr_up$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
qastr_dn <- bitr(qastr_dn$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

qcocul_up <- bitr(qcocul_up$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
qcocul_dn <- bitr(qcocul_dn$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

bg_genes <- bitr(bg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

gene_list2test <- list("Astro-Up" = qastr_up$ENTREZID, "Cocult-Up" = qcocul_up$ENTREZID, 
                       "Astro-Down" = qastr_dn$ENTREZID, "Cocult-Down" = qcocul_dn$ENTREZID)
comp_clust_reac <- compareCluster(geneCluster = gene_list2test, 
                                  universe = bg_genes$ENTREZID,
                                  pAdjustMethod = "BH",
                                  fun =  enrichPathway,
                                  organism = "human", 
                                  readable=T)


outfile_png <- paste0( outPath_reac, "/bothcult_clust_cmp_enrichPath_up_down_reg.png")

###
#saving the object
head(comp_clust_reac)
save(comp_clust_reac, file = gsub(".png$", ".Rda", outfile_png) )

#saving results
enr_dt <- as.data.table(comp_clust_reac )
ofile <- gsub(".png$", ".tsv", outfile_png)
fwrite(enr_dt , file = ofile , quote=F, sep="\t")

head(comp_clust_reac, 2)
title <- paste0("Compare cluster-enrich Reactome Pathways ")
png(outfile_png, width = 12, height = 12, units = "in", res=300)
print(dotplot(comp_clust_reac, font.size = 12, showCategory=30) + ggtitle(title) +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.25, hjust=0.25)) ) # orderBy = "pvalue" ,, decreasing=F - error with enrichPathway

dev.off()


### compareClust() astr coculture overlapping and independent
#############################################################
outPath_overlapInde <- paste0(out_res_clustProf, "dge_overlap_inde/")
if(! dir.exists(outPath_overlapInde)){
  dir.create(outPath_overlapInde)
}
#fdr <5 % & FC >=2
outFile <- paste0(outPath_overlapInde,"/dge_fdr5_fc2overlap_inde.png")
get_overlap_inde_run_upsetr_clustProf(dge_all[abs(lfcShrink_log2FoldChange) >= fc2_abs_LFC_cutoff & padj < 0.05, ],unique(dge_all$gene_name), outFile)



## Chekcing the enrichment of up and down regulated genes at FC 1.5, 2.5 and 3 - Investigating if gene ratio is more in coculture than astrocytyes
#################################################################################

#running compareCluster() for the astr and cocult up and down togehter.
#code above put in a function
run_compareCluster4astr_cocul_up_dn <- function(dge_astr_dt, dge_cocul_dt, abs_LFC_cutoff, outfile_png ){
  qastr_up <- dge_astr_dt[ padj < 0.05 &
                             abs(lfcShrink_log2FoldChange) >= abs_LFC_cutoff &
                             reg == "up", .(gene_name)]
  qastr_dn <- dge_astr_dt[ padj < 0.05 &
                             abs(lfcShrink_log2FoldChange) >= abs_LFC_cutoff &
                             reg == "down", .(gene_name) ]
  
  qcocul_up <- dge_cocul_dt[ padj < 0.05 &
                               abs(lfcShrink_log2FoldChange) >= abs_LFC_cutoff &
                               reg == "up",  .(gene_name)]
  qcocul_dn <- dge_cocul_dt[ padj < 0.05 &
                               abs(lfcShrink_log2FoldChange) >= abs_LFC_cutoff &
                               reg == "down",  .(gene_name)]
  
  bg_genes <- unique(dge_astr_dt$gene_name, dge_cocul_dt$gene_name) 
  
  qastr_up <- bitr(qastr_up$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  qastr_dn <- bitr(qastr_dn$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  qcocul_up <- bitr(qcocul_up$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  qcocul_dn <- bitr(qcocul_dn$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  bg_genes <- bitr(bg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  gene_list2test <- list("Astro-Up" = qastr_up$ENTREZID, "Cocult-Up" = qcocul_up$ENTREZID, 
                         "Astro-Down" = qastr_dn$ENTREZID, "Cocult-Down" = qcocul_dn$ENTREZID)
  comp_clust_reac <- compareCluster(geneCluster = gene_list2test, 
                                    universe = bg_genes$ENTREZID,
                                    pAdjustMethod = "BH",
                                    fun =  enrichPathway,
                                    organism = "human", 
                                    readable=T)
  
  
  
  ###
  #saving the object
  head(comp_clust_reac)
  save(comp_clust_reac, file = gsub(".png$", ".Rda", outfile_png) )
  
  #saving results
  enr_dt <- as.data.table(comp_clust_reac )
  ofile <- gsub(".png$", ".tsv", outfile_png)
  fwrite(enr_dt , file = ofile , quote=F, sep="\t")
  
  head(comp_clust_reac, 2)
  title <- paste0("Compare cluster-enrich Reactome Pathways ")
  png(outfile_png, width = 12, height = 12, units = "in", res=300)
  print(dotplot(comp_clust_reac, font.size = 12, showCategory=30) + ggtitle(title) +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.25, hjust=0.25)) ) 
  
  dev.off()
  
}

run_compareCluster4astr_cocul_up <- function(dge_astr_dt, dge_cocul_dt, abs_LFC_cutoff, outfile_png ){
  qastr_up <- dge_astr_dt[ padj < 0.05 &
                             abs(lfcShrink_log2FoldChange) >= abs_LFC_cutoff &
                             reg == "up", .(gene_name)]
  
  qcocul_up <- dge_cocul_dt[ padj < 0.05 &
                               abs(lfcShrink_log2FoldChange) >= abs_LFC_cutoff &
                               reg == "up",  .(gene_name)]
  
  bg_genes <- unique(dge_astr_dt$gene_name, dge_cocul_dt$gene_name) 
  
  qastr_up <- bitr(qastr_up$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  qcocul_up <- bitr(qcocul_up$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  bg_genes <- bitr(bg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  

  ##only up
  gene_list2test <- list("Astro-Up" = qastr_up$ENTREZID, "Cocult-Up" = qcocul_up$ENTREZID)
  comp_clust_reac <- compareCluster(geneCluster = gene_list2test, 
                                    universe = bg_genes$ENTREZID,
                                    pAdjustMethod = "BH",
                                    fun =  enrichPathway,
                                    organism = "human", 
                                    readable=T)
  
  
  
  ###
  #saving the object
  save(comp_clust_reac, file = gsub(".png$", ".Rda", outfile_png) )
  
  #saving results
  enr_dt <- as.data.table(comp_clust_reac )
  ofile <- gsub(".png$", ".tsv", outfile_png)
  fwrite(enr_dt , file = ofile , quote=F, sep="\t")
  
  title <- paste0("Compare cluster-enrich Reactome Pathways ")
  png(outfile_png, width = 12, height = 12, units = "in", res=300)
  print(dotplot(comp_clust_reac, font.size = 12, showCategory=30) + ggtitle(title) +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.25, hjust=0.25)) ) 
  
  dev.off()
  
}

#ora GO -ont

run_compareCluster_oraGo4ont_astr_cocul_up_dn <- function(dge_astr_dt, dge_cocul_dt, abs_LFC_cutoff, outfile_prefix, ont2check= c("BP", "CC") ){
  #'compare cluster run for the up & down reg of both cultures and 
  #only the up reg of both'
  
  qastr_up <- dge_astr_dt[ padj < 0.05 &
                             abs(lfcShrink_log2FoldChange) >= abs_LFC_cutoff &
                             reg == "up", .(gene_name)]
  qastr_dn <- dge_astr_dt[ padj < 0.05 &
                             abs(lfcShrink_log2FoldChange) >= abs_LFC_cutoff &
                             reg == "down", .(gene_name) ]
  
  qcocul_up <- dge_cocul_dt[ padj < 0.05 &
                               abs(lfcShrink_log2FoldChange) >= abs_LFC_cutoff &
                               reg == "up",  .(gene_name)]
  qcocul_dn <- dge_cocul_dt[ padj < 0.05 &
                               abs(lfcShrink_log2FoldChange) >= abs_LFC_cutoff &
                               reg == "down",  .(gene_name)]
  
  bg_genes <- unique(dge_astr_dt$gene_name, dge_cocul_dt$gene_name) 
  
  
  gene_list2test <- list("Astro-Up" = qastr_up$gene_name, "Cocult-Up" = qcocul_up$gene_name, 
                         "Astro-Down" = qastr_dn$gene_name, "Cocult-Down" = qcocul_dn$gene_name)
  
  gene_list2test_up <- list("Astro-Up" = qastr_up$gene_name, "Cocult-Up" = qcocul_up$gene_name )
  
  
  gene_list2test_astr <- list("Up regulated genes" = qastr_up$gene_name, "Down regulated genes" = qastr_dn$gene_name )
  
  
  gene_list2test_cocul <- list("Up regulated genes" = qcocul_up$gene_name, "Down regulated genes" = qcocul_dn$gene_name)
  
  
  #compare cluster GO 4 ontology
  ret_ont <- lapply(ont2check, function(ontology){
    o_file <- paste0(outfile_prefix, "_up_dn_reg", ontology, ".png")
    ret <- run_compareCluster_go4ontology_genelist(gene_list2test, bg_genes, ontology,  o_file, title_suffix="")
    
    o_file <- paste0(outfile_prefix, "_up_reg", ontology, ".png")
    ret <- run_compareCluster_go4ontology_genelist(gene_list2test_up, bg_genes, ontology,  o_file, title_suffix="")

    o_file <- paste0(gsub("bothcult", "astr", outfile_prefix), "_up_down_reg", ontology, ".png")
    ret <- run_compareCluster_go4ontology_genelist(gene_list2test_astr, unique(dge_astr_dt$gene_name), ontology,  o_file, title_suffix="")
        
    o_file <- paste0(gsub("bothcult", "cocul", outfile_prefix), "_up_down_reg", ontology, ".png")
    ret <- run_compareCluster_go4ontology_genelist(gene_list2test_cocul, unique(dge_cocul_dt$gene_name), ontology,  o_file, title_suffix="")
    
  })
  return(T)

}


# log2(1.5) = 0.5849625
#log2(2.5) = 1.321928
# log2(3) = 1.584963

#abs_LFC_cutoff <-c(0.6, 1, 1.3, 1.6) #FC, 1.5, 2, 2.5, 3
abs_LFC_cutoff  <- 1
dge_astr_all <- dge_all[ comparison== comparison_astr, ]
dge_cocult_all <- dge_all[ comparison== comparison_cocul, ]

bg_astr_genes <- dge_astr_all$gene_name
bg_cocult_genes <- dge_cocult_all$gene_name

ret <- lapply(abs_LFC_cutoff, function(lfc){
  #lfc <- abs_LFC_cutoff[1]
  if( lfc == 1) {
    outPath_reac <- paste0(out_res_clustProf, "dge/reactome/")
    outPath_go <- paste0(out_res_clustProf, "dge/go/")
  } else {
    outPath_reac <- paste0(out_res_clustProf, "dge_absLFC", sub("\\.", "pt", lfc) , "/reactome/")
    outPath_go <- paste0(out_res_clustProf, "dge_absLFC", sub("\\.", "pt", lfc) , "/go/")
    
  }
  
  if(!dir.exists(outPath_go)){
    dir.create(outPath_go, recursive = T)
    dir.create(outPath_reac, recursive = T)
  }
  
  #reac
  # outfile_pref <- paste0( outPath_reac,"/astr")
  # run_reac_ora_compclust_4up_dn_reg_genes( dge_astr_all, abs_LFC_cutoff, outfile_pref, cmpclust_listnameprefix = "Astrocytes")
  # 
  # outfile_pref <- paste0( outPath_reac,"/cocul")
  # run_reac_ora_compclust_4up_dn_reg_genes( dge_cocult_all, abs_LFC_cutoff, outfile_pref, cmpclust_listnameprefix = "Coculture")
  # 
  # 
  # out_png <- paste0( outPath_reac, "/bothcult_clust_cmp_enrichPath_up_down_reg.png")
  # #running compareCluster() for the astr and cocult up and down togehter.
  # run_compareCluster4astr_cocul_up_dn(dge_astr_all, dge_cocult_all, abs_LFC_cutoff, out_png )
  # 
  #comparecluster with up only
  # out_png <- paste0( outPath_reac, "/bothcult_clust_cmp_enrichPath_up_reg.png")
  # run_compareCluster4astr_cocul_up(dge_astr_all, dge_cocult_all, abs_LFC_cutoff, out_png )
 
  
  #GO
  out_file_pref <- paste0( outPath_go, "/bothcult_clust_cmp_enrichGO")
  run_compareCluster_oraGo4ont_astr_cocul_up_dn(dge_astr_all, dge_cocult_all, abs_LFC_cutoff, out_file_pref, ont2check= c("BP", "CC") )
  
})


#generate bar plots to view the gene ratios for the up regulated
#################################################################

library(enrichplot)
library(forcats)
library(here)
source(here("./R/plot_setting.R"))

#library(ggnewscale)
abs_LFC_cutoff <-c(0.6, 1, 1.3, 1.6) #FC, 1.5,2, 2.5, 3
cultures <- c("astr", "cocul")

eval_generatio <- function(x){ return(eval(parse(text=x)))}

enr_lfc <- lapply(abs_LFC_cutoff, function(lfc){
  #outPath_go <- paste0(out_res_clustProf, "dge/go/")
  #lfc <- abs_LFC_cutoff[1]
  if(lfc == 1){
    outPath_reac <- paste0(out_res_clustProf, "dge/reactome/")
  } else{
    outPath_reac <- paste0(out_res_clustProf, "dge_absLFC", sub("\\.", "pt", lfc) , "/reactome/")
  }
  enr <- lapply(cultures, function(cult){
    #cult <- cultures[1]
    rda_file <- paste0(outPath_reac, cult, "_ora_up_reg.Rda")
    load(rda_file)
  outpng <- sub(".Rda", "_barplot.png", rda_file)
    title <- "" # basename(rda_file)
    #plot for the culture
    n_show_category <- 30
    gen_save_bar_plot_generatio(e_ora_reac, title, outpng, n_show_category, plot_font_size =10)
      
    enr_dt <- as.data.table(e_ora_reac )
    enr_dt$culture <- cult
    #n_show_category <- 10
    return(head(enr_dt, n_show_category))
  })  
  enr <- rbindlist(enr)
  ##generating a single plot comparing both cultures
  head(enr, 2)
  
  enr[, GeneRatioValue := mapply(eval_generatio, GeneRatio)]
  setorder(enr, -culture, GeneRatioValue) 
  ggplot(data=enr, aes(x=fct_inorder(Description), y=GeneRatioValue, fill=culture)) +
    geom_bar(stat="identity", position=position_dodge())  +
    labs(x = "Enriched terms", y = "Gene Ratio") + coord_flip() +
    scale_fill_discrete(
      limits = c("astr", "cocul"),
      labels = c("Astrocyte", "Co-culture")
    )
    ggsave( filename = paste0(outPath_reac, "both_ora_up_reg_barplot.png"), dpi = IMAGE_DPI_COLOR ) 
  return(enr)
})
#enr_lfc <- rbindlist(enr_lfc)

rm(enr_lfc, outFile, enr_dt, gene_list2test, q_astr_genes, q_cocult_genes, qastr_dn, qastr_up, qcocul_dn, qcocul_up, original_gene_list)


#generating the gene ratio bar plot using the simplified GO results for the LFC1 as decided to stay with FC2.
outPath_go <- paste0(out_res_clustProf, "dge/go/")
n_show_category <- 30

ont2check <- "BP"
enr <- lapply(ont2check, function(ont){
  fname <-paste0(outPath_go, "bothcult_clust_cmp_enrichGO_up_reg", ont, "_simplify.tsv")
  
  enr_dt <- fread(fname)
  head(enr_dt, 2)
  #get the top categories
  cultures <- unique(enr_dt$Cluster)
  enr <- lapply(cultures, function(cult){
    return(head(enr_dt[ Cluster== cult, ], n_show_category))
  })
  enr_dt <- rbindlist(enr)
  rm(enr)  
  
  enr_dt[, GeneRatioValue := mapply(eval_generatio, GeneRatio)]
  
  common_terms <- enr_dt[, .N, by=.(Description)][ N > 1 ] #get teh terms common to both
  
  enr_dt_common <- enr_dt[ Description %in% common_terms$Description, ]
  
  setorder(enr_dt_common, Cluster, GeneRatioValue) 
  ggplot(data=enr_dt_common, aes(x=fct_inorder(Description), y=GeneRatioValue, fill=Cluster)) +
    geom_bar(stat="identity", position=position_dodge())  +
    labs(x = "Enriched terms", y = "Gene Ratio") + coord_flip() +
    scale_fill_discrete(
      limits = c("Astro-Up", "Cocult-Up"),
      labels = c("Astrocyte", "Co-culture")
    )
  ggsave( filename =  gsub(".tsv", "_barplot.png", fname), width = 8, height = 10, dpi = IMAGE_DPI_COLOR ) 
  
  
})  

#cnetplots

#cult <- cultures[1]
rda_file <- paste0(outPath_go, "bothcult_clust_cmp_enrichGO_up_reg", ont, ".Rda")
load(rda_file)

cnetplot(comp_clust_go, showCategory = 30)!!!!!! Need to work it right
######
#DS
######
#load the data
load(here("../results/leafcutter_additionalFilters_postdeconv/diff_splicing/combGrp_clust_effSize.rda"))

rm(cluster_recs, intron_eff_recs)
cluster_recs_single_gene_intron[ comparison == "astrocyte_treatment_vs_no_treatment", comparison := "astrocyte_treated_vs_untreated"]

cluster_recs_single_gene_intron[ comparison == "astrocyte_neuron_treatment_vs_no_treatment", comparison := "coculture_treated_vs_untreated"]

comp <- unique(cluster_recs_single_gene_intron$comparison)
#cluster_recs_single_gene_intron <- cluster_recs_single_gene_intron[ status == "Success" ,] #337908

cluster_recs_single_gene_intron <- cluster_recs_single_gene_intron[ status == "Success" & !(is.na(genes)),] #319199

cluster_recs_single_gene_intron[, abs_deltapsi := abs(cluster_recs_single_gene_intron$deltapsi)]

abs_dpsi_cutoff <- 0.1
ds_genes <- lapply(comp, function(cmp){
    #cmp <- comp[1]
    #head(cluster_recs_single_gene_intron, 2)  
    recs <- cluster_recs_single_gene_intron[ comparison== cmp,]
    #head(recs, 2)
    #select genes of interest - same filter as used to select query gene
    recs <- recs[ p.adjust < 0.05 & abs_deltapsi >= abs_dpsi_cutoff & !is.na(genes), ]
    min_recs <- recs[ , .SD[which.min(p.adjust)], by = "genes"] #the 
    return(min_recs)
  })
  ds_genes <- rbindlist(ds_genes)
  table(ds_genes$comp)
  
  
  out_res_clustProf_ds <- paste0(out_res_clustProf, "ds", gsub("0.", "_pt", abs_dpsi_cutoff), "/")
  outPath_ds_go <- paste0(out_res_clustProf_ds, "go/")
  outPath_ds_reac <- paste0(out_res_clustProf_ds, "reac/")
  
  if(!dir.exists(out_res_clustProf_ds)){
    dir.create(out_res_clustProf_ds)
    dir.create(outPath_ds_reac)
    dir.create(outPath_ds_go)
  }
  
  
#ORA -GO

q_ds_astr_genes <- ds_genes[ comparison== comparison_astr, genes ]
bg_ds_astr_genes <- unique(cluster_recs_single_gene_intron[comparison== comparison_astr & !is.na(genes), genes])
q_ds_cocul_genes <- ds_genes[ comparison== comparison_cocul, genes ]
bg_ds_cocul_genes <- unique(cluster_recs_single_gene_intron[comparison== comparison_cocul & !is.na(genes), genes])


#Running for all ont - giving an error. Therefore exlcuding it
outfile_pref <- paste0( outPath_ds_go,"/astr_ora")
run_go_ora_exclAll(q_ds_astr_genes , bg_ds_astr_genes, outfile_pref, n_show_category=30,  png_width =10, png_ht=10, png_units="in", plot_font_size =12, png_res=300) #runs the enrichGO() for each ontology folloewd by simplify()
# function also runs enrichGO() for 'ALL' ontology - simplify() does not work on this

outfile_pref <- paste0( outPath_ds_go,"/cocul_ora")
run_go_ora_exclAll(q_ds_cocul_genes, bg_ds_cocul_genes, outfile_pref, n_show_category=30,  png_width =10, png_ht=10, png_units="in", plot_font_size =12, png_res=300) #runs the enrichGO() for each ontology folloewd by simplify()



##ORA - reac
outfile_pref <- paste0( outPath_ds_reac,"/astr_ora")
run_reac_ora (q_ds_astr_genes, bg_ds_astr_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300, plot_title_suffix = "-DS genes in astrocytes")

outfile_pref <- paste0( outPath_ds_reac,"/cocul_ora")
run_reac_ora (q_ds_cocul_genes, bg_ds_cocul_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300, plot_title_suffix =  "-DS genes in coculture")


#COMPARE CLUSTER ASTR VS COCUL
#function(qgene_set1, qset1_listname, qgene_set2, qset2_listname, bg_genes, outfile_png, title_suffix="")
run_compareCluster_reac(q_ds_astr_genes, qset1_listname ="Astrocyte",q_ds_cocul_genes, qset2_listname="Coculture",unique( bg_ds_astr_genes, bg_ds_cocul_genes), 
                   outfile_png = paste0( outPath_ds_reac, "/bothcult_clust_cmp_enrichPath.png"), title_suffix="")

#NOTE
#For the ORA, Reactome pathhways - no enrihcmet terms in astro and only 1 in cocult. Not much seen in GO too
#doing  gsea, ordering by dpsi

#GSEA -GO   -all genes tested
head(ds_genes, 2)

min_recs_dpsi <- cluster_recs_single_gene_intron[ , .SD[which.min(deltapsi)], by = c("comparison", "genes")] #the 
#ds_astr_all <- min_recs_dpsi[ comparison== comparison_astr & ! is.na(genes), ]
#original_gene_list <- ds_astr_all[ comparison== comparison_astr , deltapsi]
#names(original_gene_list) <- ds_astr_all$genes
#run_gsea4genelist(original_gene_list,  paste0(outPath_ds_go, "astr_gse_go.png") , n_show_cate = 30 ) #Error in get("GO2ONT", envir = GO_DATA) : object 'GO2ONT' not found 



ontologies <- c("BP", "CC", "MF")
ret <- lapply(ontologies, function(o){
  #o<- ontologies[1]
  #astr
  ##go
  cat("\n ont= ", o, "\n")
  ds_astr_all <- min_recs_dpsi[ comparison== comparison_astr & ! is.na(genes), ]
  original_gene_list <- ds_astr_all[ comparison== comparison_astr , deltapsi]
  names(original_gene_list) <- ds_astr_all$genes
  l_out_png <- paste0(outPath_ds_go, "astr_gse_", o,".png") 
  run_gsea4genelist_4ont(original_gene_list, l_out_png, n_show_cate = 30, ont_name=o )
  rm(ds_astr_all, l_out_png, original_gene_list)
  
  #cocult
  ##go
  ds_cocul_all <- min_recs_dpsi[ comparison== comparison_cocul & ! is.na(genes), ]
  original_gene_list <- ds_cocul_all[ comparison== comparison_cocul , deltapsi]
  names(original_gene_list) <- ds_cocul_all$genes
  l_out_png <- paste0(outPath_ds_go, "cocul_gse_", o,".png") 
  run_gsea4genelist_4ont(original_gene_list, l_out_png, n_show_cate = 30, ont_name=o )
  rm(ds_cocul_all, l_out_png, original_gene_list)
  
})

head(min_recs_dpsi)

##reac
ds_astr_all <- min_recs_dpsi[ comparison== comparison_astr & ! is.na(genes), ]
l_astr_all <- copy(ds_astr_all)
head(l_astr_all, 2)
setnames(l_astr_all, old= c("genes", "deltapsi" ), new=c("gene_symbol", "avg_log2FC"))  
names(l_astr_all)
run_gsea_reac(l_astr_all, paste0(outPath_ds_reac, "astr_gse_reac.png") , n_show_cate = 30, png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300, nPermSimple_val= 100000 ) #as suggested by clustprof, no enrichment at default

rm(ds_astr_all)

ds_cocul_all <- min_recs_dpsi[ comparison== comparison_cocul & ! is.na(genes), ]
l_cocul_all <- copy(ds_cocul_all)
setnames(l_cocul_all, old= c("genes", "deltapsi" ), new=c("gene_symbol", "avg_log2FC"))  
names(l_cocul_all)
run_gsea_reac(l_cocul_all, paste0(outPath_ds_reac, "cocul_gse_reac.png") , n_show_cate = 30, png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300, nPermSimple_val= 100000 ) #as suggested by clustprof, no enrichment at default

rm(ds_cocul_all)



#compare DS, DGE
######################
bg_astr_genes <- dge_astr_all$gene_name
bg_cocult_genes <- dge_cocult_all$gene_name

q_astr_genes <- dge_astr_all[ padj < 0.05 &
                                abs(lfcShrink_log2FoldChange) >= fc2_abs_LFC_cutoff, gene_name]
q_cocult_genes <- dge_cocult_all[ padj < 0.05 &
                                    abs(lfcShrink_log2FoldChange) >= fc2_abs_LFC_cutoff, gene_name]


#comparecluster - enrichGO

run_compareCluster_go(q_astr_genes, "DGE", q_ds_astr_genes, "DS", bg_astr_dge_ds, outfile_png= paste0( out_res_clustProf, "dge_ds/go/astr_clust_cmp_dge_ds_enrichGo.png"), title_suffix=" - Astrocytes")

run_compareCluster_go(q_cocult_genes, "DGE", q_ds_cocul_genes, "DS", bg_cocul_dge_ds, outfile_png= paste0( out_res_clustProf, "dge_ds/go/cocul_clust_cmp_dge_ds_enrichGo.png"), title_suffix=" - Coculture")

#comparecluster - enrichPathway
run_compareCluster_reac(q_astr_genes, "DGE", q_ds_astr_genes, "DS", bg_astr_dge_ds, outfile_png= paste0( out_res_clustProf, "dge_ds/reac/astr_clust_cmp_dge_ds_enrichPath.png"), title_suffix=" - Astrocytes")

run_compareCluster_reac(q_cocult_genes, "DGE", q_ds_cocul_genes, "DS", bg_cocul_dge_ds, outfile_png= paste0( out_res_clustProf, "dge_ds/reac/cocul_clust_cmp_dge_ds_enrichPath.png"), title_suffix=" - Coculture")



##DS - running the enrichment for ORA- reactome varying the abs dpsi threshold

abs_dpsis <- 0.1 #c(0.05, 0.1, 0.15, 0.2, 0.25)# c( 0.30, 0.35)  # c(0.05, 0.15, 0.2) #0.1,
retval <- lapply( abs_dpsis, function(abs_dpsi_cutoff){
  #abs_dpsi_cutoff <- abs_dpsis[1]
    ds_genes <- lapply(comp, function(cmp){
      #cmp <- comp[1]
      #head(cluster_recs_single_gene_intron, 2)  
      recs <- cluster_recs_single_gene_intron[ comparison== cmp,]
      #head(recs, 2)
      #select genes of interest - same filter as used to select query gene
      recs <- recs[ p.adjust < 0.05 & abs_deltapsi >= abs_dpsi_cutoff & !is.na(genes), ]
      min_recs <- recs[ , .SD[which.min(p.adjust)], by = "genes"] #the 
      return(min_recs)
    })
    ds_genes <- rbindlist(ds_genes)
    table(ds_genes$comp)
    
    
    out_res_clustProf_ds <- paste0(out_res_clustProf, "ds", sub("0.", "_pt", abs_dpsi_cutoff), "/")
    outPath_ds_go <- paste0(out_res_clustProf_ds, "go/")
    outPath_ds_reac <- paste0(out_res_clustProf_ds, "reac/")
    outPath_dge_ds <- paste0(out_res_clustProf, "dge_ds", sub("0.", "_pt", abs_dpsi_cutoff), "/")
    
    if(!dir.exists(out_res_clustProf_ds)){
      dir.create(out_res_clustProf_ds)
      dir.create(outPath_ds_reac)
      dir.create(outPath_ds_go)
      dir.create(paste0(outPath_dge_ds, "reac/"), recursive = T)
      dir.create(paste0(outPath_dge_ds, "go/"), recursive = T)
      dir.create(paste0(outPath_dge_ds, "kegg/"), recursive = T)
      
    }
    
    q_ds_astr_genes <- ds_genes[ comparison== comparison_astr, genes ]
    bg_ds_astr_genes <- unique(cluster_recs_single_gene_intron[comparison== comparison_astr & !is.na(genes), genes])
    q_ds_cocul_genes <- ds_genes[ comparison== comparison_cocul, genes ]
    bg_ds_cocul_genes <- unique(cluster_recs_single_gene_intron[comparison== comparison_cocul & !is.na(genes), genes])
    

    ##ORA - reac - defaul geneset sizes
   #  outfile_pref <- paste0( outPath_ds_reac,"/astr_ora")
   #  run_reac_ora (q_ds_astr_genes, bg_ds_astr_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300, plot_title_suffix = "-DS genes in astrocytes")
   # 
   #  outfile_pref <- paste0( outPath_ds_reac,"/cocul_ora")
   # run_reac_ora (q_ds_cocul_genes, bg_ds_cocul_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300, plot_title_suffix =  "-DS genes in coculture")

      

   ##ORA - reac - min geneset size=5
  #  minGeneSetSize <- 2 #5
  #  outfile_pref <- paste0( outPath_ds_reac,"/astr_ora_minGS", minGeneSetSize)
  #  run_reac_ora (q_ds_astr_genes, bg_ds_astr_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300, plot_title_suffix = "-DS genes in astrocytes", minGeneSetSize)
  # 
  #  outfile_pref <- paste0( outPath_ds_reac,"/cocul_ora_minGS", minGeneSetSize)
  #  run_reac_ora (q_ds_cocul_genes, bg_ds_cocul_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300, plot_title_suffix =  "-DS genes in coculture", minGeneSetSize)
  # return()
  #  
   #set maxgenesize
   # maxGeneSetSize <- 1000
   # outfile_pref <- paste0( outPath_ds_reac,"/astr_ora_minGS", minGeneSetSize, "_maxGS", maxGeneSetSize )
   # run_reac_ora (q_ds_astr_genes, bg_ds_astr_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300, plot_title_suffix = "-DS genes in astrocytes", minGeneSetSize , maxGeneSetSize)
   # 
   # outfile_pref <- paste0( outPath_ds_reac,"/cocul_ora_minGS", minGeneSetSize, "_maxGS", maxGeneSetSize )
   # run_reac_ora (q_ds_cocul_genes, bg_ds_cocul_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300, plot_title_suffix =  "-DS genes in coculture", minGeneSetSize, maxGeneSetSize)
   # return()
     
    #compare dge, ds
    bg_astr_dge_ds <- unique(bg_ds_astr_genes, bg_astr_genes)
    bg_cocul_dge_ds <- unique(bg_ds_cocul_genes, bg_cocult_genes)
    
  
    #comparecluster - enrichGO
    run_compareCluster_go(q_astr_genes, "DGE", q_ds_astr_genes, "DS", bg_astr_dge_ds, outfile_png= paste0( outPath_dge_ds,"/go/astr_clust_cmp_dge_ds_enrichGO.png"), title_suffix=" - Astrocytes")

    run_compareCluster_go(q_cocult_genes, "DGE", q_ds_cocul_genes, "DS", bg_cocul_dge_ds, outfile_png= paste0( outPath_dge_ds ,"/go/cocul_clust_cmp_dge_ds_enrichGO.png"), title_suffix=" - Coculture")
  
    
    #compare cluster GO 4 ontology
    ont2check <- c("BP", "CC")
    ret_ont <- lapply(ont2check, function(ontology){
    
      run_compareCluster_go4ontology(q_astr_genes, "DGE", q_ds_astr_genes, "DS", bg_astr_dge_ds, ontology,  outfile_png= paste0( outPath_dge_ds,"/go/astr_clust_cmp_dge_ds_enrichGO_",ontology, ".png"), title_suffix=" - Astrocytes")
      # 
      run_compareCluster_go4ontology(q_cocult_genes, "DGE", q_ds_cocul_genes, "DS", bg_cocul_dge_ds, ontology, outfile_png= paste0( outPath_dge_ds ,"/go/cocul_clust_cmp_dge_ds_enrichGO_", ontology, ".png"), title_suffix=" - Coculture")
    })
    return()
    
    #run ORA - GO
    outfile_pref <- paste0( outPath_ds_go,"/astr_ora")
    run_go_ora_exclAll (q_ds_astr_genes, bg_ds_astr_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300)
    
    outfile_pref <- paste0( outPath_ds_go,"/cocul_ora")
    run_go_ora_exclAll (q_ds_cocul_genes, bg_ds_cocul_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300)
    
    
    #comparecluster - groupGO
    #is designed for gene classification based on GO distribution at a specific level
    # ont2check <- c("BP", "CC")
    # ret_ont <- lapply(ont2check, function(ont){
    #           run_compareCluster_groupgo(q_astr_genes, "DGE", q_ds_astr_genes, "DS", ont_name = ont, outfile_png= paste0( outPath_dge_ds,"/go/astr_clust_cmp_dge_ds_grpGO_", ont, ".png"), title_suffix=" - Astrocytes")
    #     
    #       run_compareCluster_groupgo(q_cocult_genes, "DGE", q_ds_cocul_genes, "DS", ont_name = ont, outfile_png= paste0( outPath_dge_ds,"/go/cocul_clust_cmp_dge_ds_grpGO_", ont, ".png"), title_suffix=" - Coculture")
    #   
    # })    
})

# DGE DS-Running comparecluster() enrichpathways - to get all terms - no pval cutoff applied here
####################

#same code as above but excludin ght commented and setting the pval cutoff to 1
#NO enrichment observed in the astr ds genes
abs_dpsis <- 0.1 #c(0.05, 0.1, 0.15, 0.2, 0.25)# c( 0.30, 0.35)  # c(0.05, 0.15, 0.2) #0.1,
pval <- 1
out_res_clustProf_pval <- paste0(out_res_clustProf, "pvalcutoff1/")
dir.create(out_res_clustProf_pval)
retval <- lapply( abs_dpsis, function(abs_dpsi_cutoff){
  #abs_dpsi_cutoff <- abs_dpsis[1]
  ds_genes <- lapply(comp, function(cmp){
    #cmp <- comp[1]
    #head(cluster_recs_single_gene_intron, 2)  
    recs <- cluster_recs_single_gene_intron[ comparison== cmp,]
    #head(recs, 2)
    #select genes of interest - same filter as used to select query gene
    recs <- recs[ p.adjust < 0.05 & abs_deltapsi >= abs_dpsi_cutoff & !is.na(genes), ]
    min_recs <- recs[ , .SD[which.min(p.adjust)], by = "genes"] #the 
    return(min_recs)
  })
  ds_genes <- rbindlist(ds_genes)
  table(ds_genes$comp)

  out_res_clustProf_ds <- paste0(out_res_clustProf_pval, "ds", sub("0.", "_pt", abs_dpsi_cutoff), "/")
  outPath_ds_go <- paste0(out_res_clustProf_ds, "go/")
  outPath_ds_reac <- paste0(out_res_clustProf_ds, "reac/")
  outPath_dge_ds <- paste0(out_res_clustProf_pval, "dge_ds", sub("0.", "_pt", abs_dpsi_cutoff), "/")
  
  if(!dir.exists(out_res_clustProf_ds)){
    dir.create(out_res_clustProf_ds)
    dir.create(outPath_ds_reac)
    dir.create(outPath_ds_go)
    dir.create(paste0(outPath_dge_ds, "reac/"), recursive = T)
    dir.create(paste0(outPath_dge_ds, "go/"), recursive = T)
    dir.create(paste0(outPath_dge_ds, "kegg/"), recursive = T)
  }
  
  q_ds_astr_genes <- ds_genes[ comparison== comparison_astr, genes ]
  bg_ds_astr_genes <- unique(cluster_recs_single_gene_intron[comparison== comparison_astr & !is.na(genes), genes])
  q_ds_cocul_genes <- ds_genes[ comparison== comparison_cocul, genes ]
  bg_ds_cocul_genes <- unique(cluster_recs_single_gene_intron[comparison== comparison_cocul & !is.na(genes), genes])
  
  
  #run ORA - GO
   outfile_pref <- paste0( outPath_ds_go,"/astr_ora")
   run_go_ora_exclAll (q_ds_astr_genes, bg_ds_astr_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300, pval)
  
    outfile_pref <- paste0( outPath_ds_go,"/cocul_ora")
   run_go_ora_exclAll (q_ds_cocul_genes, bg_ds_cocul_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300, pval)

  return()
  
  ##ORA - reac - defaul geneset sizes
  outfile_pref <- paste0( outPath_ds_reac,"/astr_ora")
  run_reac_ora (q_ds_astr_genes, bg_ds_astr_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300, plot_title_suffix = "-DS genes in astrocytes (no pval cutoff)", pval)
  # 
  outfile_pref <- paste0( outPath_ds_reac,"/cocul_ora")
  run_reac_ora (q_ds_cocul_genes, bg_ds_cocul_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300, plot_title_suffix =  "-DS genes in coculture (no pval cutoff)", pval)
  
  
  
  #compare dge, ds
  bg_astr_dge_ds <- unique(bg_ds_astr_genes, bg_astr_genes)
  bg_cocul_dge_ds <- unique(bg_ds_cocul_genes, bg_cocult_genes)
  
  
  #comparecluster - enrichPathway
  
  run_compareCluster_reac(q_astr_genes, "DGE", q_ds_astr_genes, "DS", bg_astr_dge_ds, outfile_png= paste0( outPath_dge_ds,"/reac/astr_clust_cmp_dge_ds_enrichPath.png"), title_suffix=" - Astrocytes", pval)
  
  run_compareCluster_reac(q_cocult_genes, "DGE", q_ds_cocul_genes, "DS", bg_cocul_dge_ds, outfile_png= paste0( outPath_dge_ds ,"/reac/cocul_clust_cmp_dge_ds_enrichPath.png"), title_suffix=" - Coculture", pval)
})
