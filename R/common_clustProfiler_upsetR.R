
library(clusterProfiler)

run_compareCluster_go<- function(qgene_set1, qset1_listname, qgene_set2, qset2_listname, bg_genes, outfile_png, title_suffix="", p_val_cutoff=0.05){
  gene_list2test <- list(qgene_set1,  qgene_set2)
  names(gene_list2test) <- c(qset1_listname, qset2_listname)
  
  comp_clust_go <- compareCluster(geneCluster = gene_list2test, 
                                  universe = bg_genes,
                                  OrgDb = org.Hs.eg.db,
                                  pAdjustMethod = "BH",
                                  fun =  "enrichGO",  keyType = "SYMBOL", 
                                  pvalueCutoff= p_val_cutoff 
                                  )
  
  save(comp_clust_go, file = gsub(".png$", ".Rda", outfile_png) )
  
  #saving results
  enr_dt <- as.data.table(comp_clust_go )
  ofile <- gsub(".png$", ".tsv", outfile_png)
  fwrite(enr_dt , file = ofile , quote=F, sep="\t")

  title <- paste0("Compare cluster-enrich GO ", title_suffix)
  png(outfile_png, width = 12, height = 12, units = "in", res=300)
  print(dotplot(comp_clust_go, font.size = 12, showCategory=30) + ggtitle(title) +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.25, hjust=0.25)) ) # orderBy = "pvalue" , decreasing=F - error
  
  dev.off()
  return(T)
  
}

run_compareCluster_go4ontology<- function(qgene_set1, qset1_listname, qgene_set2, qset2_listname, bg_genes, ontology, outfile_png, title_suffix="", p_val_cutoff=0.05){
  
  gene_list2test <- list(qgene_set1,  qgene_set2)
  names(gene_list2test) <- c(qset1_listname, qset2_listname)
  
  comp_clust_go <- compareCluster(geneCluster = gene_list2test, 
                                  universe = bg_genes,
                                  OrgDb = org.Hs.eg.db,
                                  pAdjustMethod = "BH",
                                  fun =  "enrichGO",  keyType = "SYMBOL", 
                                  ont= ontology, 
                                  pvalueCutoff= p_val_cutoff)
  
  save(comp_clust_go, file = gsub(".png$", ".Rda", outfile_png) )
  
  #saving results
  enr_dt <- as.data.table(comp_clust_go )
  ofile <- gsub(".png$", ".tsv", outfile_png)
  fwrite(enr_dt , file = ofile , quote=F, sep="\t")
  
  
  title <- paste0("Compare cluster-enrich GO -",ontology, ", ", title_suffix)
  png(outfile_png, width = 12, height = 12, units = "in", res=300)
  print(dotplot(comp_clust_go, font.size = 12, showCategory=30) + ggtitle(title) +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.25, hjust=0.25)) ) 
  
  dev.off()
  return(T)
  
}

run_compareCluster_go4ontology_genelist<- function(gene_list2test, bg_genes, ontology, outfile_png, title_suffix="", p_val_cutoff=0.05, n_show_cate=30){
  
  comp_clust_go <- compareCluster(geneCluster = gene_list2test, 
                                  universe = bg_genes,
                                  OrgDb = org.Hs.eg.db,
                                  pAdjustMethod = "BH",
                                  fun =  "enrichGO",  keyType = "SYMBOL", 
                                  ont= ontology, 
                                  pvalueCutoff= p_val_cutoff)
  
  save(comp_clust_go, file = gsub(".png$", ".Rda", outfile_png) )
  
  #saving results
  enr_dt <- as.data.table(comp_clust_go )
  ofile <- gsub(".png$", ".tsv", outfile_png)
  fwrite(enr_dt , file = ofile , quote=F, sep="\t")
  
  
  title <- paste0("Compare cluster-enrich GO -",ontology, ", ", title_suffix)
  png(outfile_png, width = 12, height = 12, units = "in", res=300)
  print(dotplot(comp_clust_go, font.size = 12, showCategory=n_show_cate) + ggtitle(title) +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.25, hjust=0.25)) ) 
  
  dev.off()
  
  #simplify
  simp_enrich_res <- clusterProfiler::simplify(comp_clust_go) # using defaults
  
  png( gsub(".png", "_simplify.png", outfile_png), width = 12, height = 12, units = "in", res=300, type="cairo")
  print( dotplot(simp_enrich_res, font.size = 12, showCategory=n_show_cate) + ggtitle(paste0(title, " - Simplify")) +
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
  dev.off()
  
  simp_cluster_summary <- as.data.table( simp_enrich_res )
  fwrite(simp_cluster_summary , file = gsub(".png$", "_simplify.tsv", outfile_png) , quote=F, sep="\t")
  
  return(T)
  
}


run_compareCluster_groupgo<- function(qgene_set1, qset1_listname, qgene_set2, qset2_listname,  ont_name, outfile_png, go_level = 2, title_suffix=""){
  #groupGO() #is designed for gene classification based on GO distribution at a specific level
  #Functional Profile of a gene set at specific GO level. Given a vector of genes, groupGO() will return the GO profile at a specific level.
  gene_list2test <- list(qgene_set1,  qgene_set2)
  names(gene_list2test) <- c(qset1_listname, qset2_listname)
  
  comp_clust_go <- compareCluster(geneCluster = gene_list2test, 
                                  level = go_level,
                                  OrgDb = org.Hs.eg.db,
                                  fun =  "groupGO",  
                                  ont =ont_name, 
                                  keyType = "SYMBOL")
  
  save(comp_clust_go, file = gsub(".png$", ".Rda", outfile_png) )
  
  #saving results
  enr_dt <- as.data.table(comp_clust_go )
  ofile <- gsub(".png$", ".tsv", outfile_png)
  fwrite(enr_dt , file = ofile , quote=F, sep="\t")
  
  
  title <- paste0("Compare cluster-group GO: ",ont_name,  title_suffix)
  png(outfile_png, width = 12, height = 12, units = "in", res=300)
  print(dotplot(comp_clust_go, font.size = 12, showCategory=30)+
          ggtitle(title) +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.25, hjust=0.25)) ) 
  dev.off()
  return(T)
  
}

run_compareCluster_reac<- function(qgene_set1, qset1_listname, qgene_set2, qset2_listname, bg_genes, outfile_png, title_suffix="", p_val_cutoff=0.05){
  #gene_set1 = gene symbols
  
  entr_set1 <- bitr(qgene_set1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  entr_set2 <- bitr(qgene_set2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  entr_bggenes <- bitr(bg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  gene_list2test <- list(entr_set1$ENTREZID,  entr_set2$ENTREZID)
  names(gene_list2test) <- c(qset1_listname, qset2_listname)
  comp_clust_reac <- compareCluster(geneCluster = gene_list2test, 
                                    universe = entr_bggenes$ENTREZID,
                                    pAdjustMethod = "BH",
                                    fun =  enrichPathway,
                                    organism = "human",
                                    readable=T,
                                    pvalueCutoff= p_val_cutoff)
  
  ###
  #saving the object
  save(comp_clust_reac, file = gsub(".png$", ".Rda", outfile_png) )
  
  #saving results
  enr_dt <- as.data.table(comp_clust_reac )
  ofile <- gsub(".png$", ".tsv", outfile_png)
  fwrite(enr_dt , file = ofile , quote=F, sep="\t")
  
  head(comp_clust_reac, 2)
  title <- paste0("Compare cluster-enrich Reactome Pathways ", title_suffix)
  png(outfile_png, width = 12, height = 12, units = "in", res=300)
  print(dotplot(comp_clust_reac, font.size = 12, showCategory=30) + ggtitle(title) +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.25, hjust=0.25)) )
  
  dev.off()
  return(T)
  
}

run_gsea4genelist <- function(genelist, outfile_png, n_show_cate = 30 ){
  #function runs the gseGO() for all ontologies
 
  # omit any NA values 
  gene_list<-na.omit(genelist)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  head(gene_list); head(original_gene_list)
  
  gse <- gseGO(geneList=gene_list, 
               ont ="ALL", 
               OrgDb        = org.Hs.eg.db,
               keyType = "SYMBOL", 
               minGSSize = 10,
               maxGSSize = 500,
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               pAdjustMethod = "BH", 
               eps = 0 # clustprof msgFor some pathways, in reality P-values are less than 1e-10. You can set the `eps` argument to zero for better estimation.
               #nPermSimple = 100000  # as suggested by warning message on intial run
  ) 

  
  clusterProf_simplify_genPlot_saveResults(enrich_res = gse, ofile_png = outfile_png, n_show_category=n_show_cate,  png_width =12, png_ht=12, png_units="in", plot_font_size =12, png_res=300, simplify= F ) #simplify only applied to output from gseGO and enrichGO... 
  
}

run_gsea4genelist_4ont <- function(genelist, outfile_png, n_show_cate = 30, ont_name="BP"){
  #function runs the gseGO() for all ontologies
  #for dGE gene_list can contain the following:
  #we want the log2 fold change 
  #original_gene_list <- dt_gene_lfc$avg_log2FC
  #names(original_gene_list) <- dt_gene_lfc$gene_symbol
  
  # omit any NA values 
  gene_list<-na.omit(genelist)
  head(gene_list)
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  #head(gene_list); head(original_gene_list)
  
  gse <- gseGO(geneList=gene_list, 
               ont =ont_name, 
               OrgDb        = org.Hs.eg.db,
               keyType = "SYMBOL", 
               minGSSize = 10,
               maxGSSize = 500,
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               pAdjustMethod = "BH", 
               eps = 0 # For some pathways, in reality P-values are less than 1e-10. You can set the `eps` argument to zero for better estimation.
               #nPermSimple = 100000  # as suggested by warning message on intial run
  ) 
  head(gse, 2)
  clusterProf_simplify_genPlot_saveResults(enrich_res = gse, ofile_png = outfile_png, n_show_category=n_show_cate,  png_width =12, png_ht=12, png_units="in", plot_font_size =12, png_res=300, simplify=T ) #simplify only applied to output from gseGO and enrichGO...  
  
}

run_go_ora_exclAll <- function(q_genes, bg_genes, ofile_pref, n_show_category=30,  png_width =10, png_ht=10, png_units="in", plot_font_size =8, png_res=300, p_val_cutoff =0.05){
  #ofile_pref <- outfile_pref;q_genes <- q_genes$gene_symbol
  #function run the enrichGO() for each ontology folloewd by simplify()
  # function also runs enrichGO() for 'ALL' ontology - simplify() does not work on this
  # the enrich result sobjects are saved as rda, and tsv
  #dot plots generated
  
  ofile_png <- paste0(ofile_pref, ".png")

  
  ego_mf <- enrichGO(gene          = q_genes,
                     universe      = bg_genes,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     keyType = "SYMBOL", 
                     pvalueCutoff = p_val_cutoff)
  
  simp_ego_mf <- clusterProfiler::simplify(ego_mf) # using defaults
  
  
  ego_cc <- enrichGO(gene          = q_genes,
                     universe      = bg_genes,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     keyType = "SYMBOL", 
                     pvalueCutoff = p_val_cutoff)
  simp_ego_cc <- clusterProfiler::simplify(ego_cc) # using defaults
  
  
  ego_bp <- enrichGO(gene          = q_genes,
                     universe      = bg_genes,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     keyType = "SYMBOL", 
                     pvalueCutoff = p_val_cutoff)
  simp_ego_bp <- clusterProfiler::simplify(ego_bp) # using defaults
  
  #saving the objects
  enrichgo_obj_list <- list("ego_bp"= ego_bp, "simp_ego_bp" = simp_ego_bp, 
                            "ego_cc" =ego_cc, "simp_ego_cc" = simp_ego_cc, 
                            "ego_mf" =ego_mf, "simp_ego_mf"= simp_ego_mf)
  names(enrichgo_obj_list)
  save(enrichgo_obj_list, file = gsub(".png$", ".Rda", ofile_png) )
  
  #saving simplified results
  simp_dt <- as.data.table(simp_ego_bp )
  ofile_bp <- gsub(".png$", "_enrGO_BP_simplify.tsv", ofile_png)
  fwrite(simp_dt , file = ofile_bp , quote=F, sep="\t")
  
  simp_dt <- as.data.table(simp_ego_cc )
  ofile_cc <-  gsub(".png$", "_enrGO_CC_simplify.tsv", ofile_png) 
  fwrite(simp_dt , file = ofile_cc , quote=F, sep="\t")
  
  simp_dt <- as.data.table(simp_ego_mf )
  ofile_mf <- gsub(".png$", "_enrGO_MF_simplify.tsv", ofile_png) 
  fwrite(simp_dt , file = ofile_mf, quote=F, sep="\t")
  
  if(nrow(as.data.table(ego_bp)) >0){
    res <- gen_save_dot_plot(ego_bp, title= "enrich GO - BP", 
                             gsub("_simplify.tsv$", ".png", ofile_bp) , n_show_category, 
                             png_width, png_ht , png_units , plot_font_size, png_res)
  }
  if(nrow(as.data.table(ego_cc)) >0){
    res <- gen_save_dot_plot(ego_cc, title= "enrich GO - CC", 
                             gsub("_simplify.tsv$", ".png", ofile_cc) , n_show_category, 
                             png_width, png_ht , png_units , plot_font_size, png_res)
  }
  if(nrow(as.data.table(ego_mf)) >0){
    res <- gen_save_dot_plot(ego_mf, title= "enrich GO - MF", 
                             gsub("_simplify.tsv$", ".png", ofile_mf), n_show_category, 
                             png_width, png_ht , png_units , plot_font_size, png_res)
  }

  #simplify plots
  if(nrow(as.data.table(simp_ego_bp)) >0){
      res <- gen_save_dot_plot(simp_ego_bp, title= "enrich GO - BP", 
                           gsub(".tsv$", ".png", ofile_bp) , n_show_category, 
                           png_width, png_ht , png_units , plot_font_size, png_res)
  }
  if(nrow(as.data.table(ego_cc)) >0){
    res <- gen_save_dot_plot(simp_ego_cc, title= "enrich GO - CC", 
                           gsub(".tsv$", ".png", ofile_cc) , n_show_category, 
                           png_width, png_ht , png_units , plot_font_size, png_res)
  }
  if(nrow(as.data.table(ego_mf)) >0){
    res <- gen_save_dot_plot(simp_ego_mf, title= "enrich GO - MF", 
                           gsub(".tsv$", ".png", ofile_mf), n_show_category, 
                           png_width, png_ht , png_units , plot_font_size, png_res)
  }
    return(T)
}


clusterProf_simplify_genPlot_saveResults <- function(enrich_res, ofile_png, n_show_category=20,  png_width =10, png_ht=10, png_units="in", plot_font_size =8, png_res=300, simplify=T ){
  #ofile_png <-  #outfile_png; enrich_res <- gse; simplify <- F
  #head(gse[ , 1:11])
  if(simplify){ 
    simp_enrich_res <- clusterProfiler::simplify(enrich_res) # using defaults
    
    png( gsub(".png", "_simplify.png", ofile_png), width = png_width, height = png_ht, units = png_units, res=png_res, type="cairo")
    print( dotplot(simp_enrich_res, font.size = plot_font_size, showCategory=n_show_category, orderBy = "pvalue", decreasing=F) +
             theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
    dev.off()
    
    simp_cluster_summary <- as.data.table(simp_enrich_res )
    fwrite(simp_cluster_summary , file = gsub(".png$", "_simplify.tsv", ofile_png) , quote=F, sep="\t")
    
    png(ofile_png, width = png_width, height = png_ht, units = png_units, res=png_res, type="cairo")
    print( dotplot(enrich_res, font.size = plot_font_size, showCategory=n_show_category, orderBy = "pvalue", decreasing=F) +
             theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
    dev.off()
  } else {
    png(ofile_png, width = png_width, height = png_ht, units = png_units, res=png_res, type="cairo")
    print( dotplot(enrich_res, font.size = plot_font_size, showCategory=n_show_category) +
             theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
    dev.off()
  }  
  
  ##to save the enrichments
  cluster_summary <- as.data.table(enrich_res )
  
  fwrite(cluster_summary , file = gsub(".png$", ".tsv", ofile_png) , quote=F, sep="\t")
  if(simplify){
      save(enrich_res, simp_enrich_res, file = gsub(".png$", ".Rda", ofile_png) )
  } else {
    save(enrich_res, file = gsub(".png$", ".Rda", ofile_png) )
    
  }
  return(T)
  
  
}

run_reac_ora <- function(q_genes, bg_genes, ofile_pref, n_show_category=30,  png_width =10, png_ht=10, png_units="in", plot_font_size =8, png_res=300, plot_title_suffix="", minGeneSetSize = 10, maxGeneSetSize = 500, p_val_cutoff=0.05){
  #ofile_pref <- outfile_pref;q_genes <- q_genes$gene_symbol
  #function run the enrichPathway() for each ontology folloewd by simplify()
  #  the enrich result sobjects are saved as rda, and tsv
  #dot plots generated
  
  #q_genes <- q_ds_astr_genes
  #bg_genes <- bg_ds_astr_genes
  #p_val_cutoff <- pval
  
  ofile_png <- paste0(ofile_pref, ".png")
  entr_qgenes <- bitr(q_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  entr_bggenes <- bitr(bg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  nrow(entr_qgenes); length(q_genes)
  nrow(entr_bggenes); length(bg_genes)
  e_ora_reac <- enrichPathway(gene          = entr_qgenes$ENTREZID,
                              universe      = entr_bggenes$ENTREZID,
                              organism = "human",
                              minGSSize = minGeneSetSize,
                              maxGSSize = maxGeneSetSize,
                              pAdjustMethod = "BH", pvalueCutoff = p_val_cutoff, readable = TRUE
                              )
  #head(entr_qgenes$ENTREZID, 2); head(entr_bggenes$ENTREZID, 2) 
  #getwd()
  #fwrite(entr_qgenes,"astr_ds_qgenes.tsv", quote=F, sep="\t"  )
  #saving the object
  save(e_ora_reac, file = gsub(".png$", ".Rda", ofile_png) )
  
  #saving results
  enr_dt <- as.data.table(e_ora_reac )
  ofile_bp <- gsub(".png$", "_enrREAC.tsv", ofile_png)
  fwrite(enr_dt , file = ofile_bp , quote=F, sep="\t")
  
  if(nrow(enr_dt) >0){
    res <- gen_save_dot_plot(e_ora_reac, title= paste0("enrich Reactome Pathways ",plot_title_suffix) , 
                             gsub(".tsv$", ".png", ofile_bp) , n_show_category, 
                             png_width, png_ht , png_units , plot_font_size, png_res)
  }  
  return(T)
}

run_gsea_reac <- function(dt_gene_lfc, ofile_png,n_show_category=30,  png_width =10, png_ht=10, png_units="in", plot_font_size =8, png_res=300, nPermSimple_val=1000){
  #function runs the gsePathway() 
  #dt_gene_lfc <-l_astr_all #allgenes_dt
  #nPermSimple_val set to clusterprofilers default
  
  entr_genes <- bitr(dt_gene_lfc$gene_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  head(entr_genes)
  dt_entr <- merge(dt_gene_lfc, entr_genes, by.x="gene_symbol", by.y="SYMBOL")
  #we want the log2 fold change 
  original_gene_list <- dt_entr$avg_log2FC
  names(original_gene_list) <- dt_entr$ENTREZID
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  head(gene_list); head(original_gene_list)
  
  gse_reac <-  gsePathway(geneList = gene_list, 
                          organism = "human", 
                          pAdjustMethod = "BH", 
                          verbose = FALSE,
                          eps=0, 
                          nPermSimple =nPermSimple_val)
  
  # png(outfile_png, width= 10, height = 10, units = "in", res= 300)
  #  dotplot(gse, showCategory=30, orderBy = "pvalue")
  #  dev.off()
  
  #outfile_png <- paste0( outPath,"/clustProf_cmp_enrichGO_scDGE_up_down_reg.png")
  
  #saving the object
  save(gse_reac, file = gsub(".png$", ".Rda", ofile_png) )
  
  #saving results
  enr_dt <- as.data.table(gse_reac )
  ofile_bp <- gsub(".png$", ".tsv", ofile_png)
  fwrite(enr_dt , file = ofile_bp , quote=F, sep="\t")
  
  
  res <- gen_save_dot_plot(gse_reac, title= "GSEA- Reactome Pathways", 
                           gsub(".tsv$", ".png", ofile_bp) , n_show_category, 
                           png_width, png_ht , png_units , plot_font_size, png_res)
  
  return(T)
}

gen_save_dot_plot <- function(enrich_res, title, ofile_png, n_show_category=20,  png_width =10, png_ht=10, png_units="in", plot_font_size =8, png_res=300){
  
  #enrich_res <- ego_mf; ofile_png <-  gsub("_simplify.tsv$", ".png", ofile_mf)
#  png(ofile_png, width = png_width, height = png_ht, units = png_units, res=png_res, type="cairo")
  png(ofile_png, width = png_width, height = png_ht, units = png_units, res=png_res)
  print(dotplot(enrich_res, font.size = plot_font_size, showCategory=n_show_category, orderBy = "pvalue", decreasing=F) + ggtitle(title))
  
  dev.off()
  return(ofile_png)
}

gen_save_bar_plot_generatio <- function(enrich_res, title, ofile_png, n_show_category=20,  png_width =10, png_ht=10, png_units="in", plot_font_size =8, png_res=300){
  
  #enrich_res <- ego_mf; ofile_png <-  gsub("_simplify.tsv$", ".png", ofile_mf)
  #  png(ofile_png, width = png_width, height = png_ht, units = png_units, res=png_res, type="cairo")
  png(ofile_png, width = png_width, height = png_ht, units = png_units, res=png_res)
  print(barplot(enrich_res,  x = "GeneRatio", font.size = plot_font_size, showCategory=n_show_category, orderBy = "pvalue", decreasing=F) + 
          ggtitle(title)) + 
          labs(x= "Gene Ratio", y = "Enriched terms")
  
  dev.off()
  return(ofile_png)
}


generate_save_upsetRPlot <-function(genes_list, outFname,x_label){
  p <-upset(fromList(genes_list), sets = names(genes_list), keep.order = TRUE, nintersects = 25, order.by = "freq",text.scale = c(1,1,1,1,1.5,1.5 ), 
            sets.x.label = x_label)
  png(outFname,width =8, height = 8, units="in", res=300, type= "cairo") #
  print(p)
  dev.off()
  return(T)
  
}

get_overlap_inde_run_upsetr_clustProf <- function(dge_dt, bg_genes, ofile_png){
  #dge_dt <- dge genes query set;  ofile_png <- outFileUpsetR_up; upsetr_xLabel <- "Treated vs untreated - up regulated"
  
  #function genertes the upsetr plots and runs compareCluster for enrichGO and enrichPathway -
  #for the overlapping and independendent up and down regualted genes in astr and cocult DGE
  
  #get the overlaps and unique genes
  astr_up <- dge_dt[ (comparison== "astrocyte_treated_vs_untreated" & reg== "up") ,.(gene_name, lfcShrink_log2FoldChange)]
  cocul_up <- dge_dt[( comparison== "coculture_treated_vs_untreated" & reg== "up") ,.(gene_name, lfcShrink_log2FoldChange)]
  
  astr_dn <- dge_dt[ (comparison== "astrocyte_treated_vs_untreated" & reg== "down") ,.(gene_name, lfcShrink_log2FoldChange)]
  cocul_dn <- dge_dt[( comparison== "coculture_treated_vs_untreated" & reg== "down") ,.(gene_name, lfcShrink_log2FoldChange)]
  
  upgenes <- list("Astrocyte" = astr_up$gene_name, "Co-culture" = cocul_up$gene_name )
  dngenes <- list("Astrocyte" = astr_dn$gene_name, "Co-culture" = cocul_dn$gene_name )
  
  #generate upsetR plot
  generate_save_upsetRPlot( upgenes, gsub(".png", "_up_upsetr.png", ofile_png),x_label="Treated vs untreated - up regulated" )
  generate_save_upsetRPlot( dngenes, gsub(".png", "_down_upsetr.png", ofile_png),x_label="Treated vs untreated - down regulated" )
  
  overlap_up <- intersect(astr_up$gene_name, cocul_up$gene_name)
  overlap_dn <- intersect(astr_dn$gene_name, cocul_dn$gene_name) 
  
  unique_astr_up <-setdiff(astr_up$gene_name, cocul_up$gene_name)
  unique_astr_dn <-setdiff(astr_dn$gene_name, cocul_dn$gene_name)
  unique_cocult_up <-setdiff(cocul_up$gene_name, astr_up$gene_name)
  unique_cocult_dn <-setdiff(cocul_dn$gene_name, astr_dn$gene_name)
  
  gene_list2test <- list("Astrocyte independent-up" = unique_astr_up, 
                         "Co-culture independent-up" = unique_cocult_up,
                         "Astrocyte independent-down" = unique_astr_dn,
                         "Co-culture independent-down" = unique_cocult_dn,
                         "Overlapping-up" = overlap_up,
                         "Overlapping-down" = overlap_dn)
  comp_clust_go <- compareCluster(geneCluster = gene_list2test, 
                               universe = bg_genes,
                               OrgDb = org.Hs.eg.db,
                               pAdjustMethod = "BH",
                               fun =  "enrichGO",  keyType = "SYMBOL")
  
  clusprof_png <- gsub( ".png", "_up_dn_clustProf_cmp_enrichGO.png", ofile_png) 
  clusterProf_simplify_genPlot_saveResults(comp_clust_go, ofile_png = clusprof_png, n_show_category=20,  png_width =12, png_ht=12, png_units="in", plot_font_size =10, png_res=300, simplify = F )
  rm(comp_clust_go)
  
  #for ontologies BP, CC & MF
  ont2check= c("BP", "CC", "MF") 
  ret_ont <- lapply(ont2check, function(ontology){
    o_file <- gsub(".png", paste0("_up_dn_clustProf_cmp_GO_", ontology, ".png"), ofile_png)
    ret <- run_compareCluster_go4ontology_genelist(gene_list2test, bg_genes, ontology,  o_file, title_suffix="")
    
  })
    
  #reactome Pathways
  unique_astr_up <- bitr(unique_astr_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  unique_astr_dn <- bitr(unique_astr_dn, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  unique_cocult_up <- bitr(unique_cocult_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  unique_cocult_dn <- bitr(unique_cocult_dn, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  overlap_up <- bitr(overlap_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  overlap_dn <- bitr(overlap_dn, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  bg_genes <- bitr(bg_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  gene_list2test <- list("Astrocyte independent-up" = unique_astr_up$ENTREZID, 
                         "Co-culture independent-up" = unique_cocult_up$ENTREZID,
                         "Astrocyte independent-down" = unique_astr_dn$ENTREZID,
                         "Co-culture independent-down" = unique_cocult_dn$ENTREZID,
                         "Overlapping-up" = overlap_up$ENTREZID,
                         "Overlapping-down" = overlap_dn$ENTREZID)
  comp_clust_reac <- compareCluster(geneCluster = gene_list2test, 
                                    universe = bg_genes$ENTREZID,
                                    pAdjustMethod = "BH",
                                    fun =  enrichPathway,
                                    organism = "human")
  
  clusprof_png <- gsub( ".png", "_up_dn_clustProf_cmp_enrichPathway.png", ofile_png) 
  clusterProf_simplify_genPlot_saveResults(comp_clust_reac, ofile_png = clusprof_png, n_show_category=20,  png_width =12, png_ht=12, png_units="in", plot_font_size =10, png_res=300, simplify = F )
  
  rm(enrich_dge, gene_list2test )
  return(T)
  
}

##
