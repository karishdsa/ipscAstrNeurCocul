#clusterprofiler run on the DEGs cluster0 vs 4 - sc data
library(data.table)
library(here)
library(openxlsx)

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

library(DT)

library(ggplot2)
source(here("./R/plot_setting.R"))

source(here("./R/common_clustProfiler_upsetR.R"))
source(here("./R/common_sc_dge.R"))

options(bitmapType='cairo')   

inPath <- here("../results/dge_sc/")
outPath <- paste0(inPath, "clust_prof/")
outPath_go <- paste0(outPath, "go/")
outPath_reac <- paste0(outPath, "reactome/")
outPath_wiki <- paste0(outPath, "wiki/")

if(!dir.exists(outPath_go)){
  dir.create(outPath_go, recursive = T)
  dir.create(outPath_reac)
  dir.create(outPath_wiki)
}

dge_infile <- paste0(inPath, "DGE_findMarker_lfcThresh0_cluster0_vs_4.xlsx")
#get the dge data
allgenes_dt <- fetch_sc_dge(dge_infile)

#GO enrichment analysis
######################
#run GSEA
run_gsea(allgenes_dt, paste0(outPath_go, "gse_go.png") , n_show_cate = 30 )

#GO over representation analysis
head(allgenes_dt, 2)
q_genes <- allgenes_dt[sig_fdr5_lfcgt1 == "Y", .(gene_symbol)]
bg_genes <- allgenes_dt$gene_symbol

outfile_pref <- paste0( outPath_go,"/ora")

run_go_ora(q_genes$gene_symbol, bg_genes, outfile_pref, n_show_category=30,  png_width =10, png_ht=10, png_units="in", plot_font_size =12, png_res=300) #runs the enrichGO() for each ontology folloewd by simplify()
# function also runs enrichGO() for 'ALL' ontology - simplify() does not work on this
  

#comparecluster GOfor the BP and CC ont
#################################
clus0_up <- allgenes_dt[ sig_fdr5_lfcgt1 == "Y" & reg == "up" ,.(gene_symbol)]

clus0_down <- allgenes_dt[ sig_fdr5_lfcgt1 == "Y" & reg == "down" ,.(gene_symbol)] #will be up regulated in clust4
bg_genes <- allgenes_dt$gene_symbol

gene_list2test <- list("Cluster 0 -up" = clus0_up$gene_symbol, 
                       "Custer 4 - up" = clus0_down$gene_symbol)

outfile_pref <- paste0( outPath_go,"/clust_cmp_")

ont2check <- c("BP", "CC")
ret_ont <- lapply(ont2check, function(ontology){
  o_file <- paste0(outfile_pref, "up_dn_reg_", ontology, ".png")
  ret <- run_compareCluster_go4ontology_genelist(gene_list2test, bg_genes, ontology,  o_file, title_suffix="")
  })



#Running the enrich Reactome
############################
outfile_pref <- paste0( outPath_reac,"/ora")

run_reac_ora (q_genes$gene_symbol, bg_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300)

#GSEA - reac
run_gsea_reac(allgenes_dt, paste0(outPath_reac, "gse_reac.png") , n_show_cate = 30, png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300 )

bg_genes <- allgenes_dt$gene_symbol

#Running ora - Reactome for the up and down regulated genes
#up
q_genes <- allgenes_dt[sig_fdr5_lfcgt1 == "Y" & reg == "up", .(gene_symbol)]
outfile_pref <- paste0( outPath_reac,"/ora_up_reg")
run_reac_ora (q_genes$gene_symbol, bg_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300)


#down
q_genes <- allgenes_dt[sig_fdr5_lfcgt1 == "Y" & reg == "down", .(gene_symbol)]
outfile_pref <- paste0( outPath_reac,"/ora_dn_reg")

run_reac_ora (q_genes$gene_symbol, bg_genes, outfile_pref, n_show_category=30,  png_width =14, png_ht=14, png_units="in", plot_font_size =12, png_res=300)

#running compareCluster()

clus0_up <- allgenes_dt[ sig_fdr5_lfcgt1 == "Y" & reg == "up" ,.(gene_symbol)]
entr_clus0_up <- bitr(clus0_up$gene_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

clus0_down <- allgenes_dt[ sig_fdr5_lfcgt1 == "Y" & reg == "down" ,.(gene_symbol)] #will be up regulated in clust4
entr_clus0_down <- bitr(clus0_down$gene_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

entr_bggenes <- bitr(allgenes_dt$gene_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

gene_list2test <- list("Cluster 0 -up" = entr_clus0_up$ENTREZID, 
                       "Custer 4 - up" = entr_clus0_down$ENTREZID)
comp_clust_reac <- compareCluster(geneCluster = gene_list2test, 
                              universe = entr_bggenes$ENTREZID,
                             pAdjustMethod = "BH",
                             fun =  enrichPathway,
                             organism = "human", readable=T)
outfile_png <- paste0( outPath_reac,"/clust_cmp_enrichPath_up_down_reg.png")


###
#saving the object
save(comp_clust_reac, file = gsub(".png$", ".Rda", outfile_png) )

#saving results
enr_dt <- as.data.table(comp_clust_reac )
ofile <- gsub(".png$", ".tsv", outfile_png)
fwrite(enr_dt , file = ofile , quote=F, sep="\t")

head(comp_clust_reac, 2)
title <- "Compare cluster-enrich Reactome Pathways"
png(outfile_png, width = 12, height = 12, units = "in", res=300)
print(dotplot(comp_clust_reac, font.size = 12, showCategory=30) + ggtitle(title)) # orderBy = "pvalue" ,, decreasing=F - error with enrichPathway

dev.off()


## Running ORA - GO terms for the up and down regulated genes
#############################
bg_genes <- allgenes_dt$gene_symbol

#up reg
clus0_up <- allgenes_dt[ sig_fdr5_lfcgt1 == "Y" & reg == "up" ,.(gene_symbol)]
outfile_pref <- paste0( outPath_go,"/ora_up_reg")
run_go_ora(clus0_up$gene_symbol, bg_genes, outfile_pref, n_show_category=30,  png_width =10, png_ht=10, png_units="in", plot_font_size =12, png_res=300) #runs the enrichGO() for each ontology folloewd by simplify()

#down
clus0_down <- allgenes_dt[ sig_fdr5_lfcgt1 == "Y" & reg == "down" ,.(gene_symbol)] #will be up regulated in clust4
outfile_pref <- paste0( outPath_go,"/ora_dn_reg")
run_go_ora(clus0_down$gene_symbol, bg_genes, outfile_pref, n_show_category=30,  png_width =10, png_ht=10, png_units="in", plot_font_size =12, png_res=300) #runs the enrichGO() for each ontology folloewd by simplify()


##running compareCluster()

bg_genes <- allgenes_dt$gene_symbol

up_genes <- allgenes_dt[sig_fdr5_lfcgt1 == "Y" & reg == "up", .(gene_symbol)]
dn_genes <- allgenes_dt[sig_fdr5_lfcgt1 == "Y" & reg == "down", .(gene_symbol)]


gene_list2test <- list("Cluster 0 -up" = up_genes$gene_symbol, 
                       "Custer 4 - up" = dn_genes$gene_symbol)
comp_clust_go <- compareCluster(geneCluster = gene_list2test, 
                                  universe = bg_genes,
                                  pAdjustMethod = "BH",
                                  OrgDb = org.Hs.eg.db,
                                  keyType = "SYMBOL",
                                  fun =  enrichGO)
outfile_png <- paste0( outPath_go,"/clust_cmp_enrichGO_up_down_reg.png")


###
#saving the object
save(comp_clust_go, file = gsub(".png$", ".Rda", outfile_png) )

#saving results
enr_dt <- as.data.table(comp_clust_go )
ofile <- gsub(".png$", ".tsv", outfile_png)
fwrite(enr_dt , file = ofile , quote=F, sep="\t")

head(comp_clust_go, 2)
title <- "Compare cluster-enrich GO"
png(outfile_png, width = 12, height = 12, units = "in", res=300)
print(dotplot(comp_clust_go, font.size = 12, showCategory=30) + ggtitle(title)) 

dev.off()

