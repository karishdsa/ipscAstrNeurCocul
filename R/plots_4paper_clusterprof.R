
library(data.table)
library(here)
library(org.Hs.eg.db)
library(ReactomePA)
library(clusterProfiler)
library(UpSetR)
source(here("./R/common_dge.R"))
source(here("./R/common_clustProfiler_upsetR.R"))

library(ggplot2)
source(here("./R/plot_setting_paper.R"))

options(bitmapType='cairo')   

#SC data
########
inPath <- here("../results/dge_sc/clust_prof/")
outPath <- here("../results/paper/")

if(!dir.exists(outPath)){
  dir.create(outPath, recursive = T)
}
rda_file <- paste0(inPath, "go/clust_cmp_up_dn_reg_BP.Rda")

simpl_rda <- paste0(dirname(rda_file), "/", gsub(".Rda", "_simplify.Rda", basename(rda_file)))

SIZE_AXIS_TEXT <-11
out_pdf <-  paste0(outPath, gsub(".Rda", "_scData_10cate.pdf", basename(simpl_rda)))
generate_dot_plot_pdf (simpl_rda, out_pdf,  n_show_category=10, title = "",  p_width =10, p_ht=8, plot_font_size =14, sc_data=T ) # only pdf needed - chnagin size for illust


out_png <-  paste0(outPath, gsub(".Rda", "_scData_10cate_landsc.png", basename(simpl_rda)))
generate_dot_plot_landscape(simpl_rda, out_png,  n_show_category=10, title = "",  png_width =11, png_ht=7, png_units="in", png_res=300, sc_data=T )

#pdf - from generate_dot_plot_landscape()
out_pdf <-  paste0(outPath, gsub(".Rda", "_scData_10cate_landsc.pdf", basename(simpl_rda)))
rda_file <- simpl_rda; n_show_category <- 10
load(rda_file)
#head(as.data.table(simp_enrich_res), 2)
levels(simp_enrich_res@compareClusterResult$Cluster)[levels(simp_enrich_res@compareClusterResult$Cluster) == "Cluster 0 -up"] <- "Astrocyte cluster 1"
levels(simp_enrich_res@compareClusterResult$Cluster)[levels(simp_enrich_res@compareClusterResult$Cluster) == "Custer 4 - up"] <- "Astrocyte cluster 2"

p <-  dotplot(simp_enrich_res, showCategory=n_show_category) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  coord_flip() +
  theme(legend.key.size = unit(1, "cm")) +
  theme(legend.position = "bottom", legend.box="vertical")
pdf( out_pdf, width = 12, height = 8)
print(  p )           
dev.off()
#without legend
out_pdf <-  paste0(outPath, gsub(".Rda", "_scData_10cate_landsc_noLegend.pdf", basename(simpl_rda)))
p <-  dotplot(simp_enrich_res, showCategory=n_show_category) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  coord_flip() +
  theme(legend.position = "none")
pdf( out_pdf, width = 7, height = 5.5)
print(  p )           
dev.off()



#BULK DATA
###########
out_res_clustProf <- here( "../results/cluster_profiler/" )
out_dge_path <- here("../results/dge_postdeconv/")

outPath_go <- here("../results/paper/")

if(!dir.exists(outPath_go)){
  dir.create(outPath_go, recursive = T)
}

inpath_enrich_rda <- here("../results/cluster_profiler/dge/go/")
#DGE

#GO enrichment analysis
#########################

#function to load the rda, simplify and generate the plot
run_simplify <- function(rda_file, out_rda ){
  #the rda should contain the comp_clust_go object
  load(rda_file)
  simp_enrich_res <- clusterProfiler::simplify(comp_clust_go) # using defaults
  save(simp_enrich_res, file = out_rda) 
  return(T)
  
}

generate_dot_plot_prePapertheme <- function(rda_file, outfile_png,  n_show_category=30, title = "",  png_width =10, png_ht=10, png_units="in", plot_font_size =12, png_res=300, sc_data=F ){
  #the rda should contain the comp_clust_go object
  #rda_file <- simpl_rda; outfile_png <- out_png
  load(rda_file)
  #head(as.data.table(simp_enrich_res), 2)
  if(sc_data){
    levels(simp_enrich_res@compareClusterResult$Cluster)[levels(simp_enrich_res@compareClusterResult$Cluster) == "Cluster 0 -up"] <- "Astrocyte cluster 1"
    levels(simp_enrich_res@compareClusterResult$Cluster)[levels(simp_enrich_res@compareClusterResult$Cluster) == "Custer 4 - up"] <- "Astrocyte cluster 2"
  }  
  
  png( outfile_png, width = png_width, height = png_ht, units = png_units, res=png_res, type="cairo")
  print( dotplot(simp_enrich_res, font.size = plot_font_size, showCategory=n_show_category) + ggtitle(title) +
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 legend.text = element_text(size = plot_font_size-4), 
                 legend.title  = element_text(size = plot_font_size-4)) )
  dev.off()
  pdf( gsub("png$", "pdf", outfile_png),  width = png_width*1.25, height = png_ht*1.25)
  print( dotplot(simp_enrich_res, font.size = plot_font_size, showCategory=n_show_category) + ggtitle(title) +
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                 legend.text = element_text(size = plot_font_size-4), 
                 legend.title  = element_text(size = plot_font_size-4)) )
  dev.off()
  
  simp_cluster_summary <- as.data.table( simp_enrich_res )
  fwrite(simp_cluster_summary , file = gsub(".png$", ".tsv", outfile_png) , quote=F, sep="\t")
  
  return(T)
  
}

generate_dot_plot <- function(rda_file, outfile_png,  n_show_category=30, title = "",  png_width =10, png_ht=10, png_units="in", plot_font_size =12, png_res=300, sc_data=F ){
  #the rda should contain the comp_clust_go object
  #rda_file <- simpl_rda; outfile_png <- out_png
  load(rda_file)
  if(sc_data){
    levels(simp_enrich_res@compareClusterResult$Cluster)[levels(simp_enrich_res@compareClusterResult$Cluster) == "Cluster 0 -up"] <- "Astrocyte cluster 1"
    levels(simp_enrich_res@compareClusterResult$Cluster)[levels(simp_enrich_res@compareClusterResult$Cluster) == "Custer 4 - up"] <- "Astrocyte cluster 2"
  }  
  p <- dotplot(simp_enrich_res, showCategory=n_show_category) + ggtitle(title) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1) )
  
  png( outfile_png, width = png_width, height = png_ht, units = png_units, res=png_res, type="cairo")
  print( p )           
  dev.off()
  pdf( gsub("png$", "pdf", outfile_png),  width = png_width*1.25, height = png_ht*1.25)
  print( p )
  dev.off()
  
  simp_cluster_summary <- as.data.table( simp_enrich_res )
  fwrite(simp_cluster_summary , file = gsub(".png$", ".tsv", outfile_png) , quote=F, sep="\t")
  
  return(T)
  
}

generate_dot_plot_pdf <- function(rda_file, outfile_pdf,  n_show_category=30, title = "",  p_width =10, p_ht=10, plot_font_size =12, sc_data=F, wrap_terms=F ){
  #the rda should contain the comp_clust_go object
  #rda_file <- simpl_rda; outfile_png <- out_png
  load(rda_file)
  #head(as.data.table(simp_enrich_res), 2)
  if(sc_data){
    levels(simp_enrich_res@compareClusterResult$Cluster)[levels(simp_enrich_res@compareClusterResult$Cluster) == "Cluster 0 -up"] <- "Astrocyte cluster 1"
    levels(simp_enrich_res@compareClusterResult$Cluster)[levels(simp_enrich_res@compareClusterResult$Cluster) == "Custer 4 - up"] <- "Astrocyte cluster 2"
  }  
  p <- dotplot(simp_enrich_res, showCategory=n_show_category) + ggtitle(title) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  pdf( outfile_pdf,  width = p_width, height = p_ht)
  print( p )
  dev.off()
  return(T)
}

generate_dot_plot_landscape <- function(rda_file, outfile_png,  n_show_category=30, title = "",  png_width =10, png_ht=10, png_units="in", plot_font_size =12, png_res=300, sc_data=F ){
  #the rda should contain the comp_clust_go object
  #rda_file <- simpl_rda; outfile_png <- out_png;   n_show_category <- 10
  load(rda_file)
  #head(as.data.table(simp_enrich_res), 2)
  if(sc_data){
    levels(simp_enrich_res@compareClusterResult$Cluster)[levels(simp_enrich_res@compareClusterResult$Cluster) == "Cluster 0 -up"] <- "Astrocyte cluster 1"
    levels(simp_enrich_res@compareClusterResult$Cluster)[levels(simp_enrich_res@compareClusterResult$Cluster) == "Custer 4 - up"] <- "Astrocyte cluster 2"
  }  
  
  p <-  dotplot(simp_enrich_res, showCategory=n_show_category) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    coord_flip()
  
  png( outfile_png, width = png_width, height = png_ht, units = png_units, res=png_res, type="cairo")
  print(  p )           
  dev.off()
  pdf( gsub("png$", "pdf", outfile_png),  width = png_width*1.25, height = png_ht*1.25)
  print( p)
  dev.off()
  
  simp_cluster_summary <- as.data.table( simp_enrich_res )
  fwrite(simp_cluster_summary , file = gsub(".png$", ".tsv", outfile_png) , quote=F, sep="\t")
  
  return(T)
  
}

# gnerateing plots for DGE up and down regulaed - astr and cocult
#runnign simplify() fiollowed bt a plot
rda_file <- paste0(inpath_enrich_rda, "astr_clust_cmp_enrichGO_up_down_regBP.Rda" )
simpl_rda <- paste0(inpath_enrich_rda, gsub(".Rda", "_simplify.Rda", basename(rda_file)))

#run_simplify(rda_file, simpl_rda) #already run once

out_png <-  paste0(outPath_go, gsub(".Rda", "_10cate.png", basename(simpl_rda)))
generate_dot_plot (simpl_rda, out_png,  n_show_category=10, title = "",  png_width =10, png_ht=10, png_units="in", plot_font_size =16, png_res=300 )

out_pdf <-  paste0(outPath_go, gsub(".Rda", "_10cate.pdf", basename(simpl_rda)))
SIZE_AXIS_TEXT <-11
generate_dot_plot_pdf (simpl_rda, out_pdf,  n_show_category=10, title = "",  p_width =8, p_ht=8, plot_font_size =14, sc_data=F ) # only pdf needed - chnagin size for illust
  

#for the paper - increase font size and no title - pre new theme
################################################
#
#runnign simplify() fiollowed bt a plot

rda_file <- paste0(inpath_enrich_rda, "cocul_clust_cmp_enrichGO_up_down_regBP.Rda" )
simpl_rda <- paste0(inpath_enrich_rda, gsub(".Rda", "_simplify.Rda", basename(rda_file)))
#run_simplify(rda_file, simpl_rda) #already run once

out_png <-  paste0(outPath_go, gsub(".Rda", "_10cate.png", basename(simpl_rda)))
generate_dot_plot (simpl_rda, out_png,  n_show_category=10, title = "",  png_width =9, png_ht=9, png_units="in", png_res=300 )

out_pdf <-  paste0(outPath_go, gsub(".Rda", "_10cate.pdf", basename(simpl_rda)))
SIZE_AXIS_TEXT <-11
generate_dot_plot_pdf (simpl_rda, out_pdf,  n_show_category=10, title = "",  p_width =8, p_ht=8, plot_font_size =14, sc_data=F ) # only pdf needed - chnagin size for illust

#generate bar plots to view the gene ratios for the up regulated
#################################################################

library(enrichplot)
library(forcats)
library(here)

cultures <- c("astr", "cocul")
eval_generatio <- function(x){ return(eval(parse(text=x)))}

#generating the gene ratio bar plot using the simplified GO results for the LFC1 as decided to stay with FC2.
n_show_category <- 30  #ncate 10 width =10 height=6    # ncate 30  width 10, height 8
img_width <- 5
img_height <- 6
ont2check <- "BP"

enr <- lapply(ont2check, function(ont){
  #ont <- ont2check[1]
  fname <-paste0(inpath_enrich_rda, "bothcult_clust_cmp_enrichGO_up_reg", ont, "_simplify.tsv")
  enr_dt <- fread(fname)
  head(enr_dt, 2)
  #get the top categories
  cultures <- unique(enr_dt$Cluster)
  enr <- lapply(cultures, function(cult){
    return(head(enr_dt[ Cluster== cult, ], n_show_category))
  })
  enr_dt <- rbindlist(enr)
  rm(enr)  
  #head(enr_dt, 1)
  enr_dt[, GeneRatioValue := mapply(eval_generatio, GeneRatio)]
  
  common_terms <- enr_dt[, .N, by=.(Description)][ N > 1 ] #get teh terms common to both
  
  enr_dt_common <- enr_dt[ Description %in% common_terms$Description, ]
  setnames(enr_dt_common, old= "Cluster", new = "Culture")
  
  setorder(enr_dt_common, Culture, GeneRatioValue) 
  head(enr_dt_common, 2)
  table(enr_dt_common$Culture)
  enr_dt_common[ Culture =="Astro-Up", Culture:="astrocyte_treated"] 
  enr_dt_common[ Culture =="Cocult-Up", Culture:="coculture_treated"] 
  ggplot(data=enr_dt_common, aes(x=fct_inorder(Description), y=GeneRatioValue, fill=Culture)) +
    geom_bar( stat="identity", position=position_dodge())  +
    labs(x = "Enriched terms", y = "Gene Ratio") + coord_flip() +
    scale_fill_manual(    values = palette_culture,
                          limits = c("astrocyte_treated", "coculture_treated"),
                          labels = c("Astrocyte", "Co-culture")
    ) +
    theme(   axis.line = element_line(colour = "black"))

  out_png <-  paste0(outPath_go, gsub(".tsv", "_barplot_colScheme.png", basename(fname)))
  
  ggsave( filename =  out_png, width = img_width, height = img_height, dpi = IMAGE_DPI_COLOR ) 
  ggsave( filename =  gsub("png$", "pdf", out_png), width = img_width+2, height = img_height+4 ) 
  
  
})  



# editing and DGE
#####################################

inpath_enrich_rda <- here("../results/mRNA_editing/dge_ds_overlap/cluster_profiler/check4paper/dge/go/" )

#runnign simplify() fiollowed bt a plot
rda_file <- paste0(inpath_enrich_rda, "bothCult_clust_cmp_bgSigDGE_sigin_bothDgeDed_FCcutoff_enrichGO_BP.Rda" )
simpl_rda <- paste0(inpath_enrich_rda, gsub(".Rda", "_simplify.Rda", basename(rda_file)))

#run_simplify(rda_file, simpl_rda) #already run once

outPath <- here("../results/paper/")

out_pdf <-  paste0(outPath, gsub(".Rda", "_10cate.pdf", basename(simpl_rda)))
SIZE_AXIS_TEXT <-11
generate_dot_plot_pdf (simpl_rda, out_pdf,  n_show_category=10, title = "",  p_width =8, p_ht=8, plot_font_size =14, sc_data=F ) # only pdf needed - chnagin size for illust

out_pdf <-  paste0(outPath, gsub(".Rda", "_30cate.pdf", basename(simpl_rda)))
generate_dot_plot_pdf (simpl_rda, out_pdf,  n_show_category=30, title = "",  p_width =10, p_ht=8, plot_font_size =14, sc_data=F ) # only pdf needed - chnagin size for illust


