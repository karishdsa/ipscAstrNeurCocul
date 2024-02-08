
source(here("R/common_sample.R"))
source(here("R/plot_setting.R"))

melt_dt_celltypeprop <- function( celltype_dt){
  #melts the scaden predictions/sc cell type prop dt
  mdata <- melt(celltype_dt, id=c(names(celltype_dt)[ grep("celltype_", names(celltype_dt), invert = T)] ))
  setnames(mdata, old = c("variable", "value"), new = c("cell_type", "cell_type_prop"))
  mdata$cell_type <- (gsub("celltype_", "", mdata$cell_type))  
  return(mdata)
}


get_sc_celltype_prop <- function(sc_data_path, suffix="_celltypes.txt") {
  #Function returns the cell type proportions
  #calculated using the celltypes.txt file, generated for scaden, which contains the celltype for all the cells in the sample
  #
  #sc_path <- here("../results/deconv_scaden/input/")
  sc_files <- list.files(sc_data_path, pattern = paste0("*", suffix), full.names = T)
  ret_ct_prop <- lapply(sc_files, function(f){
    #f <- sc_files[1]
    sample <- gsub(suffix, "", basename(f)) 
    l_dt <- fread(f)
    celltype_prop <- l_dt %>% 
      dplyr::group_by(Celltype) %>% 
      dplyr::summarise(n = n()) %>%
      dplyr::mutate(sample = sample,
                    scRNA_proportions = n/sum(n))  %>% 
      dplyr::select(-n) %>% 
      tidyr::spread(key = Celltype, value = scRNA_proportions)
    return(setDT(celltype_prop))     
    
  })
  ct_prop <- rbindlist(ret_ct_prop, fill=T)
  ct_prop <- ct_prop %>% 
    mutate( cell_treatment = gsub("ipsc_", "", sample)) %>%
    separate( cell_treatment, c("cells", "treatment"), sep= "_") %>%
    mutate(cells = ifelse(cells == 'astro', "astrocyte", ifelse(cells == "neur", "neuron", "co-culture"))) %>%
    mutate(treatment = ifelse(treatment == "ctrl", "untreated", "oligomer_treated"))
  
  return(ct_prop )
  
}

get_sc_celltype_prop_with_indi_info<- function(sc_celltypeprop_path){
  #returns the sc cellt ype prop along with the details of the controls used in the sc RNAseq

  celltype_prop <- get_sc_celltype_prop(sc_data_path = sc_celltypeprop_path)
  celltype_prop$var_cutoff <- "scRNA-seq"
  celltype_prop$nCells <- "scRNA-seq"
  celltype_prop$run <- "scRNA-seq"
  celltype_prop$indi <- "ctrl_1_JOM_inhouse"
  celltype_prop[ cells=="neuron", indi :="ctrl_NTN"]
  return(celltype_prop)
}

get_scaden_celltype_prop <- function(scaden_pred_file){
  #scaden_pred_file <- here("../results/deconv_scaden/predictions/scaden_predictions.txt") 
  cell_type_prop <- fread(scaden_pred_file)
  #add the sample info
  #unique(cell_type_prop$V1)
  sample_info <- get_sample_info()
  sample_info[ , V1:= paste0( cells, treatment, "_", indi) ]
  head(sample_info,2)
  sample_info <- unique(sample_info[, .(cells, indi, treatment, V1)]  )
  sample_info[  ( cells == "astrocyte_neuron")  , cells:="co-culture"]
  sample_info[  ( treatment == "no_treatment")  , treatment:="untreated"]
  sample_info[  ( treatment == "1uM_syn_oligomer")  , treatment:="oligomer_treated"]
  
  
  cell_type_prop <-  merge(cell_type_prop, sample_info, by = "V1", all.x=T)
  setnames(cell_type_prop, old="V1", new= "sample")
  return(cell_type_prop)
}

get_scaden_predictions <- function(prediction_path, pred_file_prefix, extract_nSamples_from_fname=F){
  #for output files eg scaden_pred_comb_var0_1_cells100.txt - scaden run at different var_cutoff/ ncells 
  #pred_file_prefix <- "scaden_pred_comb_" # use to search the pred files and also replaced in the file name to get the var, nCells, nSamples
  #extract_nSamples_from_fname = if the filename includes the nSamples run, split the fname further to get the nSamples in addition to var and nCells
  #calls get_scaden_celltype_prop() for all the pred files 
  #returns a data table with the cell type proportions, indi, treatment, cells, var_cutoff, nCells
  
  
  pred_files <- list.files(prediction_path, pattern=paste0(pred_file_prefix, "*"), full.names = T)
  celltype_prop <- lapply(pred_files, function(f){
    #f<- pred_files[1]
    prop <- get_scaden_celltype_prop(f)
    prop$run <- gsub( paste0(pred_file_prefix, "|\\.txt"), "",basename(f))
    return(prop)
  })
  celltype_prop <- rbindlist(celltype_prop, use.names = T)
  #head(celltype_prop,2)
  celltype_prop[, c("var_cutoff", "nCells") := tstrsplit(run, "_cells", fixed=TRUE)]
  celltype_prop$var_cutoff = gsub("_", ".", gsub("var", "", celltype_prop$var_cutoff))
  if(extract_nSamples_from_fname) {
    celltype_prop[, c("nCells", "nSamples") := tstrsplit(nCells, "_samples", fixed=TRUE)]
  }
  return(celltype_prop)
}

plot_celltype_prop <- function(ctype_dt, outPath, title){
  #ctype_dt <- rbindlist(list(sc_celltype_prop, scaden_celltype_prop), fill=T)
  
  mdata <- melt(ctype_dt, id=c(names(ctype_dt)[ grep("celltype_", names(ctype_dt), invert = T)] ))
  
  setnames(mdata, old = c("variable", "value"), new = c("cell_type", "cell_type_prop"))
  mdata[ , cell_treatment:= paste0( cells, "\n", gsub("_", "\n", treatment) ) ]
  mdata$cell_type <- (gsub("celltype_", "", mdata$cell_type))  
  
  head(mdata)
  
  theme_plot <- theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  ggplot( mdata, aes(x = cell_type, y = cell_type_prop*100, group= cell_type) ) +
    geom_boxplot( aes(col = cell_type)  ) +
    facet_grid(cells~data) +
    ggtitle(title) +
    theme_plot
  
  ggsave(paste0(outPath, "/celltype_prop_byCells.png"), dpi = 300)
  
  ggplot( mdata, aes(x = cell_type, y = cell_type_prop*100, group= cell_type) ) +
    geom_boxplot( aes(col = cell_type)  ) +
    facet_grid(cell_treatment~data) +
    ggtitle(title) +
    theme_plot
  
  ggsave(paste0(outPath, "/celltype_prop_byCellsTreatment.png"), dpi = 300)
  
  
}

getCorrelation <- function(x, y, corr_method="spearman") {
  #x <- recs$fdr_cyt; y<- recs$fdr_nuc
  #return spearman correlation by default else pearson
  #x <- m_bulk$cell_type_prop; y <- m_sc$cell_type_prop
  
  dtCorr <- data.table(method=character(), corr=numeric(), corrPval=numeric())
  if(corr_method== "spearman"){
    corr <- cor.test(x, y, method="spearman")
  } else {
    corr <- cor.test(x, y )
    
  }
  dtCorr <- rbind(dtCorr, list(corr_method, corr$estimate, formatC(corr$p.value, format = "e", digits = 2)))
  
  return(dtCorr)
}

get_corr_sc_scaden <- function( sc_dt, scaden_dt) {
  #returns the spearman correlation of the sc and scaden cell type proportions
  # for each celltye_treatment in the sc RNAseq data,  ge the samples from the bulk decon,
  #   for each var and nCell comb,
  #        get the correlation for each indi ( sample) in thebulk with that of sc
  
  
  cell_treat <- unique(sc_dt$cell_treat)
  
  corr_sc_bulk <- lapply(cell_treat, function(ct){
    # for each celltye_treatment in the sc RNAseq data,  ge the samples from the bulk decon,
    #   for each var and nCell comb,
    #        get the correlation for each indi ( sample) in thebulk with that of sc
    
    #ct <- cell_treat[1]
    l_sc <- sc_dt[ cell_treat == ct, ]
    m_sc<- melt( l_sc, id=c(names(l_sc)[ grep("celltype_", names(l_sc), invert = T)] ))
    setnames(m_sc, old = c("variable", "value"), new = c("cell_type", "cell_type_prop"))
    m_sc$cell_type <- (gsub("celltype_", "", m_sc$cell_type))  
    
    
    l_bulk <- scaden_dt[ cell_treat == ct, ]
    var_cutoffs_run <- unique(l_bulk$var_cutoff)
    nCells_run <- unique(l_bulk$nCells)
    
    ret_var_ncells <- lapply(var_cutoffs_run, function(v){
      ret_ncells <-lapply(nCells_run, function(cell_count){
        #v <- var_cutoffs_run[1]; cell_count <- nCells_run[1]  
        bulk_recs <- l_bulk[ var_cutoff ==v & nCells ==cell_count, ]
        
        ctrls <- unique(bulk_recs$indi) 
        ret_ctrl <- lapply(ctrls, function(control){
          #control <- ctrls[1]
          ctrl_bulk_recs <- bulk_recs[ indi == control, ]
          m_bulk<- melt( ctrl_bulk_recs, id=c(names(ctrl_bulk_recs)[ grep("celltype_", names(ctrl_bulk_recs), invert = T)] ))
          setnames(m_bulk, old = c("variable", "value"), new = c("cell_type", "cell_type_prop"))
          m_bulk$cell_type <- (gsub("celltype_", "", m_bulk$cell_type))  
          
          setorder(m_bulk, cell_type) 
          setorder(m_sc, cell_type) #ascending
          corr_dt <- getCorrelation(m_sc$cell_type_prop, m_bulk$cell_type_prop)
          
          corr_dt$indi <- control
          #corr_dt$scaden_run <- unique(m_bulk$run)
          return(corr_dt)
          
        })
        ret_ctrl <- rbindlist(ret_ctrl)
        ret_ctrl$ncells <- cell_count
        ret_ctrl$var_cutoff <- v
        return( ret_ctrl )     
      })
      return(rbindlist( ret_ncells) )
    })
    ret_var_ncells <- rbindlist( ret_var_ncells)
    ret_var_ncells$cell_treat <- ct
    return(ret_var_ncells)
  } )
  
  corr_sc_bulk <- rbindlist(corr_sc_bulk) 
  return(corr_sc_bulk)
  
  
}

generate_plot_ctProp_scVsScad <- function( pred_path, pred_file_prefix, sc_ctProp_path, var2sel ,  title, outFname, matchedCtrls=NULL, nCells2sel=NULL, extract_nSamples_from_fname =F, facetbyCells =T) {
  # function generates a plot for the matched controls 
  # celltype prop form the scRNA seq vs scaden for a var cut off specified
  # final if else condition to be updated if all controls to be plotted
  #if facetbyCells ==F will facet by nSamples
  
  scaden_celltype_prop <- get_scaden_predictions(pred_path, pred_file_prefix, extract_nSamples_from_fname)

  scaden_celltype_prop$run <- NULL
  scaden_celltype_prop <- scaden_celltype_prop[ var_cutoff == var2sel , ]
  if(is.null(matchedCtrls)) {
    scaden_indiFilt <- scaden_celltype_prop 
    } else {
    scaden_indiFilt <- scaden_celltype_prop[ indi %in% matchedCtrls , ] #& ncells > nCells_gt & var_sel %in% var_sel, ]
    
    }
  if(! is.null(nCells2sel)){
    scaden_indiFilt <- scaden_indiFilt[ nCells == nCells2sel, ]
  }
  

  #head(scaden_indiFilt, 2)
  scaden_mdata <- melt_dt_celltypeprop(scaden_indiFilt)
  
  #head(scaden_mdata)
  setnames(scaden_mdata , old= "cell_type_prop", new= "scaden_cellTypeProp")
  scaden_mdata$sample <- NULL
  rm(scaden_indiFilt)
  
  sc_celltype_prop <- get_sc_celltype_prop_with_indi_info(sc_celltypeprop_path = sc_ctProp_path )
  sc_mdata <- melt_dt_celltypeprop(sc_celltype_prop)
  setnames(sc_mdata , old= "cell_type_prop", new= "sc_cellTypeProp")
  sc_mdata <- unique(sc_mdata[ , .(cells, treatment, indi, cell_type, sc_cellTypeProp)])
  #head(scaden_mdata, 2); head(sc_mdata, 2)
  
  
  comb <- merge(scaden_mdata, sc_mdata, by = c("cells", "treatment", "indi", "cell_type"), all.x=T)
  head(comb, 2)
  #comb$cell_type <- paste0("Celltype ", comb$cell_type)
  #comb$nCells <- as.numeric(comb$nCells)
  head(comb$scaden_cellTypeProp)
  comb[is.na(sc_cellTypeProp), sc_cellTypeProp:=0] 
  
  length(comb$scaden_cellTypeProp)
  
  if( is.null(matchedCtrls)) {
    stop("Code to plot currently only for the matched controls. All controls will make it too busy. Code will need to be updated.")
    
  } 
  if(facetbyCells) {
    ggplot( comb, aes(x = scaden_cellTypeProp*100, y = sc_cellTypeProp*100)) +
      geom_point(aes(col=cell_type, shape = treatment)) +
      facet_grid(cells ~ as.numeric(nCells)) +
      geom_abline(slope=1, intercept = 0, linetype="dashed") +
      ggtitle(title)
    
  } else {
    ggplot( comb, aes(x = scaden_cellTypeProp*100, y = sc_cellTypeProp*100)) +
      geom_point(aes(col=cell_type, shape = treatment)) +
      facet_grid(cells ~ as.numeric(nSamples)) +
      geom_abline(slope=1, intercept = 0, linetype="dashed") +
      ggtitle(title)
    
  }
  ggsave( outFname, dpi = 300, height = 8, width = 10)
  
}

