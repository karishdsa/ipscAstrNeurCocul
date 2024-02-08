library(data.table)

#get the significant genes - functions returns the gene id's
get_sig_genes_dge <- function(deseq_rda_path, comparisons_ref_vs_base){
  sig_genes <- lapply(comparisons_ref_vs_base, function(ref_vs_base){
    #n<-1
    infile <- paste0(deseq_rda_path, "deseq2_", ref_vs_base, ".rda"  )
    
    cat("\nReading file for : ", ref_vs_base, "\n\t", infile )
    load(infile)
    
    return(sig_res$geneid)
  })
  names(sig_genes) <- comparisons_ref_vs_base
  return(sig_genes)
}


get_sig_upReg_genes_dge <- function(deseq_rda_path, comparisons_ref_vs_base){
  sig_genes <- lapply(comparisons_ref_vs_base, function(ref_vs_base){
    #n<-1
    infile <- paste0(deseq_rda_path, "deseq2_", ref_vs_base, ".rda"  )
    
    cat("\nReading file for : ", ref_vs_base, "\n\t", infile )
    load(infile)
    up <- sig_res[ which( sig_res$lfcShrink_log2FoldChange > 0), ]
    
    return(up$geneid)
  })
  names(sig_genes) <- comparisons_ref_vs_base
  
  return(sig_genes)
}


get_sig_downReg_genes_dge <- function(deseq_rda_path, comparisons_ref_vs_base){
  sig_genes <- lapply(comparisons_ref_vs_base, function(ref_vs_base){
    #n<-1
    infile <- paste0(deseq_rda_path, "deseq2_", ref_vs_base, ".rda"  )
    
    cat("\nReading file for : ", ref_vs_base, "\n\t", infile )
    load(infile)
    down <- sig_res[ which( sig_res$lfcShrink_log2FoldChange < 0), ]
    
    return(down$geneid)
  })
  names(sig_genes) <- comparisons_ref_vs_base
  
  return(sig_genes)
}


get_all_genes_dge <- function(deseq_rda_path, comparisons_ref_vs_base){
  #function returns all the result genes a a col of the sig genes are up or down regulated
  all_genes <- lapply(comparisons_ref_vs_base, function(ref_vs_base){
    #n<-1
    infile <- paste0(deseq_rda_path, "deseq2_", ref_vs_base, ".rda"  )
    cat("\nReading file for : ", ref_vs_base, "\n\t", infile )
    load(infile)
    df <- all_res[which(!is.na(all_res$padj)),] 
    df$reg <- NA
    df[which(df$padj < 0.05 & df$lfcShrink_log2FoldChange < 0), "reg"] <- "down"
    df[which(df$padj < 0.05 & df$lfcShrink_log2FoldChange > 0), "reg"] <- "up"
  
    return(df)
    })
  return(rbindlist(all_genes))
}

getGeneNameBiotype <- function(ensgIDs, chrName=FALSE) {
  #get the gene names, biotypes for the genes
  # using the file generated for the ENST ids as it also contains gene level info
  #ensgIDs <- head(rownames( vsd))
  #ensgIDs <- genes2Plot
  #note - geneInfo <- getGeneNameBiotype(genes2plot, chrName=T)
  #geneInfo <- geneInfo[ chromosome_name %in% c(1:22, "X", "Y", "MT"), ]
  
  txInfo <- fread("/home/kdsa/ipscAstrocyteNeuron/results/salmon_quant/salmon_txiQuantification_annEns93.tsv")
  if(chrName){
    return( unique( txInfo[ ensembl_gene_id %in% ensgIDs, .(ensembl_gene_id, external_gene_name, gene_biotype, chromosome_name) ]) )
  } else {
  return( unique( txInfo[ ensembl_gene_id %in% ensgIDs, .(ensembl_gene_id, external_gene_name, gene_biotype) ]) )
  }
}

getGeneIDBiotype <- function(geneNames) {
  #get the gene ids, biotypes for the genes
  # using the file generated for the ENST ids as it also contains gene level info
  #ensgIDs <- head(rownames( vsd))
  #note - geneInfo <- getGeneIDBiotype(genes2plot)
  #geneInfo <- geneInfo[ chromosome_name %in% c(1:22, "X", "Y", "MT"), ]
  
  txInfo <- fread("/home/kdsa/ipscAstrocyteNeuron/results/salmon_quant/salmon_txiQuantification_annEns93.tsv")
  return( unique( txInfo[ external_gene_name %in% geneNames, .(ensembl_gene_id, external_gene_name, gene_biotype, chromosome_name) ]) )
  
  
}

library(EnhancedVolcano)
generateVolcanoPlot<- function(df, pvalCutoff=10e-6, fcCutoff=2, outFname, title_comparison){
  #setting the default cut offs to match that set by of the package
  subTitle <- paste0(title_comparison, "\nnGenes = ", nrow(df[ which(df$pvalue< pvalCutoff & df$lfcShrink_log2FoldChange > abs(fcCutoff)),]))
  mainTitle <- paste0("log2Fold change > |", fcCutoff, "|, p-value cut-off ", pvalCutoff)
  cat(mainTitle)
  png(outFname)
  print( EnhancedVolcano(df,
                         lab = df$external_gene_name,
                         x = 'lfcShrink_log2FoldChange',
                         title= mainTitle ,
                         subtitle = subTitle,
                         y = 'pvalue',
                         FCcutoff=fcCutoff,
                         xlim = c(-12, 12),
                         pCutoff = pvalCutoff,
                         col=c('black', 'black', 'black', 'red3'), 
                         subtitleLabSize = 10)
  )
  dev.off()
  return(outFname)
}

run_generateVolcanoPlot <- function(rdaPath, comp) {
  outVolFiles <- lapply(comp, function(cs){
    print(cs)
    infile <- paste0(rdaPath, "/deseq2_", cs, ".rda")
    
    print(infile)
    load(infile)
    head(all_res, 3)
    genename <- getGeneNameBiotype(all_res$geneid)
    nrow(all_res)
    head(genename, 2)
    all_res <- merge(all_res, genename, by.x = "geneid", by.y= "ensembl_gene_id", all.x=T)
    #outFname <- paste0("/home/kdsa/ipscAstrocyteNeuron/results/dge/correct_PC23/volcano_plots/", "volPlot_", cs, ".png")
    outFname <- paste0(rdaPath, "/volcano_plots/", "volPlot_", cs, ".png")
    
    #outFname
    
    #(fold change,FC=2, Log2FC=1) 
    #in the volcano plot, the fold change indicates lfc, so passing 1 as the
    generateVolcanoPlot(all_res, pvalCutoff=10e-6, fcCutoff=1, outFname, title_comparison = cs)
    
    
    
    return(outFname)
  })
  return(unlist(outVolFiles))
}

#generate the counts using the shrunken LFC
generateDegCounts_shrunkenlfc <- function(deseq_rda_path, comparisons_ref_vs_base, outFile_tsv){
  
  cts <- lapply(comparisons_ref_vs_base, function(ref_vs_base){
    #n<-1
    infile <- paste0(deseq_rda_path, "deseq2_", ref_vs_base, ".rda"  )
    load(infile)
    
    df_counts <- data.frame(nGenes = nrow(all_res[ which(!is.na(all_res$padj)), ]), 
                            nDeg_fdr5 =nrow(sig_res), 
                            nDown_regulated = nrow( sig_res[which( sig_res$lfcShrink_log2FoldChange < 0), ]), 
                            nUp_regulated =  nrow( sig_res[which( sig_res$lfcShrink_log2FoldChange > 0), ]), 
                            comp_vs_base = ref_vs_base, 
                            stringsAsFactors=F)
    return(df_counts)
    
  })
  cts <- rbindlist(cts)
  fwrite(cts, file= outFile_tsv, sep = "\t", quote=F, row.names = F)
  
  return(cts)
}
