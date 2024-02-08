#see plot_deconv_bulk_celltype_prop.R for cell type plots with paper col scheme

library(here)
source(here("R","common_deconv.R"))
source(here("R/plot_setting_paper.R"))


outPath <- here("../results/paper/")
if(!dir.exists(outPath)){
  dir.create(outPath)
}

pred_file_prefix <- "scaden_pred_comb_"
pred_path <- here("../results/deconv_scaden/allsamples_run/predictions_nsamples/")
scaden_celltype_prop <- get_scaden_predictions(pred_path, pred_file_prefix)

#nSamples 2000 , nCells 3000, variance default 0.1
run2select <- "var0_1_cells3000_samples2000"

unique(scaden_celltype_prop$run)
#head(scaden_celltype_prop, 1)
scaden_celltype_prop <- scaden_celltype_prop[ run == run2select, ]
scaden_celltype_prop$analysis <- "Scaden"

sc_prop_path <- here("../results/deconv_scaden/input/")
sc_celltype_prop <- get_sc_celltype_prop_with_indi_info(sc_celltypeprop_path = sc_prop_path )
sc_celltype_prop$analysis <- "scRNA-seq"  

comb_dt <- rbindlist(list(sc_celltype_prop, scaden_celltype_prop), use.names = T)

head(comb_dt, 2)
setnames( comb_dt, old = c("celltype_0", "celltype_1", "celltype_2", "celltype_3", "celltype_4", "celltype_5","celltype_6", "celltype_7"), 
              new =c("Astrocyte cluster 1", "Neuron cluster 1", "Neuron cluster 2", "Neuron cluster 3", "Astrocyte cluster 2",  "Neuron cluster 4", "Neuron cluster 5", "Neuron cluster 6"))

head(comb_dt, 1)
mdata <- melt(comb_dt, id=c(names(comb_dt)[ grep("cluster", names(comb_dt), invert = T)] ))

setnames(mdata, old = c("variable", "value"), new = c("cell_type", "cell_type_prop"))
SIZE_AXIS_TEXT <- 11
SIZE_AXIS_TITLE <- SIZE_AXIS_TEXT

scale_x_discrete_cell_types <- 
  scale_x_discrete(limits = c("Astrocyte cluster 1", "Astrocyte cluster 2", "Neuron cluster 1", "Neuron cluster 2", "Neuron cluster 3",  "Neuron cluster 4", "Neuron cluster 5", "Neuron cluster 6"))

outfile_prefix <- paste0(outPath, "deconv_sc_ctype_prop_", run2select)                    
#New facet label names 
cells_labs <- c("Astrocyte", "Co-culture", "Neuron")
names(cells_labs) <- c("astrocyte", "co-culture", "neuron")


#generate the plot for only scaden deconv and coculture
#######################################################

sub_mdata <- mdata[ cells == "co-culture", ]
levels(mdata$cell_type)
theme_facet_panel <- theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background  = element_blank(),
  strip.background = element_rect(fill = "#d3d3d3", colour = "black", size = 0.5) #light grey 
)

ggplot( sub_mdata, aes(x = cell_type, y = cell_type_prop*100, fill = cell_type)) +
  geom_boxplot(fatten=1, lwd=0.5, outlier.shape =  21) + #alpha= 0.8)+
  scale_fill_manual(values = palette_clusters_comb) +
  facet_grid(cells~analysis,  labeller = labeller(cells = cells_labs)) +
  labs(x = "Cluster", y="Cell type proportion*100") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text = element_text(size = SIZE_AXIS_TEXT),
        axis.title = element_text(size = SIZE_AXIS_TITLE),
        strip.text = element_text(size = 14),
        axis.line = element_line(colour = "black")) +
  guides( fill = "none") +
  theme_facet_panel +
  scale_x_discrete_cell_types

ggsave(paste0(outfile_prefix, "_onlyCocult.png"), dpi = 300, height = 8, width = 10)
ggsave(paste0(outfile_prefix, "_onlyCocult.pdf"), height = 8, width = 11)


#only astrocyte
sub_mdata <- mdata[ cells == "astrocyte", ]
head(sub_mdata)
ggplot( sub_mdata, aes(x = cell_type, y = cell_type_prop*100, fill = cell_type)) +
  geom_boxplot(fatten=1, lwd=0.5, outlier.shape =  21) + #alpha= 0.8)+
  scale_fill_manual(values = palette_clusters_comb) +
  facet_grid(cells~analysis,  labeller = labeller(cells = cells_labs)) +
  labs(x = "Cluster", y="Cell type proportion*100") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text = element_text(size = SIZE_AXIS_TEXT),
        axis.title = element_text(size = SIZE_AXIS_TITLE),
        strip.text = element_text(size = 14),
        axis.line = element_line(colour = "black")) +
  guides(colour="none", fill = "none") +
  theme_facet_panel +
  scale_x_discrete_cell_types

ggsave(paste0(outfile_prefix, "_onlyAstr.png"), dpi = 300, height = 8, width = 10)

ggsave(paste0(outfile_prefix, "_onlyAstr.pdf"), dpi = 300, height = 8, width = 11)
