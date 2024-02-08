
# Generating the plot of the cell type proportions in the bulk, post deconvolution- for MJFF final prs, paper

#To see if there is a significant difference in the cell type proportions ( eg for celltype 0) in the treated vs untreated in each sample. 
#Analysis - deconv_bulk_celltype_prop.Rmd
#code from the notebook - updated, as including the clusternames too

library(here)
library(data.table)
library(ggplot2)

source("./R/common_sample.R")
source("./R/plot_setting_paper.R")


# load the data
celltype_propfile <-here("../results/deconv_scaden/allsamples_run/predictions_nsamples/scaden_pred_comb_var0_1_cells3000_samples2000.txt") 
outPath <- here("../results/paper/")

ct_prop <- fread(celltype_propfile)
setnames(ct_prop, old ="V1", new="sample")

sInfo <- get_sample_info()
sInfo$sample <- paste0(sInfo$cells, sInfo$treatment, "_", sInfo$indi)
#head(ct_prop, 3); head(sInfo, 2);length(unique(sInfo$sample))
sInfo <- unique( sInfo[, .(sample, cells, treatment, indi)] )
ct_prop <- merge(ct_prop, sInfo, by = "sample", all.x=T)
ct_prop[ treatment== "1uM_syn_oligomer", treatment := "treated"]
ct_prop[ treatment== "no_treatment", treatment := "untreated"]
ct_prop[ cells== "astrocyte_neuron", cells := "Co-culture"]
ct_prop[ cells== "astrocyte", cells := "Astrocyte"]

unique(ct_prop$cells)
rm(sInfo)

## Plots
# plot the proportions

mdata <- melt(ct_prop, id=c(names(ct_prop)[ grep("celltype_", names(ct_prop), invert = T)] ))
setnames(mdata, old = c("variable", "value"), new = c("cell_type", "cell_type_prop"))

#head(mdata, 3)
#excluding the neuronal samples and the cell types 5, 6 and 7 as  proportion < 0.01 was excluded in the analysis
mdata <- mdata[ cells != "neuron",]

ct2exclude<- paste0("celltype_", seq(5,7))

mdata <- mdata[ !cell_type %in% ct2exclude, ]
unique(mdata$cell_type)
#head(mdata, 3)

levels(mdata$cell_type)[levels(mdata$cell_type) == "celltype_0"] <- "Astrocyte cluster 1"
levels(mdata$cell_type)[levels(mdata$cell_type) == "celltype_4"] <- "Astrocyte cluster 2"
levels(mdata$cell_type)[levels(mdata$cell_type) == "celltype_1"] <- "Neuron cluster 1"
levels(mdata$cell_type)[levels(mdata$cell_type) == "celltype_2"] <- "Neuron cluster 2"
levels(mdata$cell_type)[levels(mdata$cell_type) == "celltype_3"] <- "Neuron cluster 3"
mdata <- droplevels(mdata)
mdata$cell_type <- factor(mdata$cell_type, levels = c("Astrocyte cluster 1", "Astrocyte cluster 2", "Neuron cluster 1", "Neuron cluster 2", "Neuron cluster 3")) #reorder

mdata$treatment <- factor(mdata$treatment, levels = c('untreated', 'treated'))
unique(mdata$treatment)


#Paper -  scheme
############################
mdata$celltype_treat <- paste0(mdata$cell_type, "_", mdata$treatment )

head(mdata, 2)
#mdata$celltype_treat <- as.factor(mdata$celltype_treat)
mdata$celltype_treat <- factor(mdata$celltype_treat , levels=c("Astrocyte cluster 1_untreated","Astrocyte cluster 1_treated",
                                                               "Astrocyte cluster 2_untreated", "Astrocyte cluster 2_treated",
                                                               "Neuron cluster 1_untreated", "Neuron cluster 1_treated",
                                                               "Neuron cluster 2_untreated", "Neuron cluster 2_treated",
                                                               "Neuron cluster 3_untreated", "Neuron cluster 3_treated"))

levels(mdata$celltype_treat)

culture <- unique(mdata$cells)
lapply(culture, function(c){
   ggplot(mdata[ cells== c,  ], aes(cell_type, cell_type_prop, fill=celltype_treat), dodge=treatment) +
    geom_boxplot(fatten=1, lwd=0.5, outlier.shape =  21) + #alpha= 0.8)+
    scale_fill_manual(values=palette_clusters) +
    labs(x= "Cell type", y= "Estimated cell type proportion") +
    theme(legend.position="") +
    theme( axis.line = element_line(colour = "black"), 
          axis.text.x = element_text( angle = 45 ,hjust = 1), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background  = element_blank()) #+
  ggsave(paste0(outPath, "deconv_bulk_prop_treatVsUntreat_", c,"_colScheme.pdf"), width =4, height = 5, units = "in", dpi =IMAGE_DPI_COLOR)
})


 ggplot( mdata, aes(x = cell_type, y = cell_type_prop, fill=celltype_treat), dodge=treatment) +
   geom_boxplot(fatten=1, lwd=0.5, outlier.shape =  21) + #alpha= 0.8)+
   scale_fill_manual(values=palette_clusters) +
   labs(x= "Cell type", y= "Estimated cell type proportion") +
   facet_wrap(~cells) +
   theme(legend.position="") +
   theme( axis.line = element_line(colour = "black"), 
          axis.text.x = element_text( angle = 45 ,hjust = 1), 
          axis.text = element_text(size = SIZE_AXIS_TEXT),
          axis.title = element_text(size = SIZE_AXIS_TITLE),
          strip.text = element_text(size = 14),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background  = element_blank(),
          strip.background = element_rect(fill = "#d3d3d3", colour = "black", size = 0.5) )#light grey ) 
 #
 ggsave(paste0(outPath, "deconv_bulk_prop_treatVsUntreat.png"), width = 8, height = 6, units = "in", dpi = IMAGE_DPI_COLOR)
 ggsave(paste0(outPath, "deconv_bulk_prop_treatVsUntreat.pdf"), width = 10, height = 8, units = "in")



