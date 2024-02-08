library(openxlsx)
fetch_sc_dge <- function(dge_xlsx_file){
  allgenes_dt <- as.data.table(read.xlsx( dge_xlsx_file))
  setnames(allgenes_dt, old="", new="gene_symbol")
  
  abs_lfc_cutoff <- 1  #great than double expression
  setDT(allgenes_dt)
  
  #allgenes_dt$sig_fdr5 <- "N"
  allgenes_dt[ p_val_adj < 0.05 , sig_fdr5:="Y" ]
  
  allgenes_dt[ p_val_adj < 0.05 & abs(avg_log2FC) > abs_lfc_cutoff, sig_fdr5_lfcgt1:="Y" ]
  allgenes_dt[sig_fdr5_lfcgt1=="Y" & avg_log2FC >0, reg := "up" ]
  allgenes_dt[sig_fdr5_lfcgt1=="Y" & avg_log2FC < 0, reg := "down" ]
  
  allgenes_dt[ p_val_adj < 0.01 , sig_fdr1:="Y" ]
  allgenes_dt[ p_val_adj < 0.01 & abs(avg_log2FC) > abs_lfc_cutoff, sig_fdr1_lfcgt1:="Y" ]
  return(allgenes_dt)
}
