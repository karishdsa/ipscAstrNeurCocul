#Visualising transcripts using ggtranscript R pacakge - https://dzhang32.github.io/ggtranscript/

# https://academic.oup.com/bioinformatics/article/38/15/3844/6617821?login=true

#Following the steps as in https://dzhang32.github.io/ggtranscript/articles/ggtranscript.html


#devtools::install_github("dzhang32/ggtranscript")
#Sys.unsetenv("GITHUB_PAT") #if err on running above line. 

library(ggtranscript)
library(rtracklayer)

library(cowplot)

library(here)
source((here("./R/plot_setting_paper.R")))
source(here("./R/common_sample.R"))
library(tidyverse)

outPath <- here("../results/paper/")

ens_ver <- 93 #106 #"93" #ensembl version
gene_of_interest <- "ADAR"  #goi =gene of interest

gtf_path <- paste0("/data/references/ensembl/gtf_gff3/v", ens_ver, "/Homo_sapiens.GRCh38.", ens_ver, ".gtf")


#gtf

#Generating the input data for ggtranscript plot from a GTF

gtf_rtrack <- rtracklayer::import(gtf_path)
gtf_rtrack <- gtf_rtrack %>% keepSeqlevels(c(1:22, "X", "Y", "MT"), pruning.mode = "coarse")

#converting to a tibble
gtf <- gtf_rtrack %>% dplyr::as_tibble()

#Filter gtf for gene of interest
goi_anno_from_gtf <-  gtf %>% 
  dplyr::filter(
    !is.na(gene_name), 
    gene_name == gene_of_interest
  ) 

# extract the required annotation columns
goi_anno_from_gtf <- goi_anno_from_gtf %>% 
  dplyr::select(
    seqnames,
    start,
    end,
    strand,
    type,
    gene_name,
    transcript_name,
    transcript_biotype
  )


#extract exons
exons <- goi_anno_from_gtf %>% dplyr::filter(type == "exon" & transcript_biotype == "protein_coding")
#head(exons, 2)

title <- "" #paste0(gene_of_interest, " -Transcripts (Ensembl v", ens_ver, ")")
plot_gene_tx <- exons %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = transcript_name
    )) +
    geom_range(
      aes(fill = transcript_biotype)
    ) +
    geom_intron(
      data = to_intron(exons, "transcript_name"),
      aes(strand = strand)
    ) +
  labs(y = "Transcript name",
       x = "")  +
  ggtitle(title) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


#Differentiating UTRs from coding seq- protein coding transcripts
#title <- "Differentiating UTRs from coding sequences in protein coding transcripts"
exons_prot_cod <- exons %>%
  dplyr::filter(transcript_biotype == "protein_coding")

# obtain cds
cds <- goi_anno_from_gtf %>% dplyr::filter(type == "CDS")

plot_utr_ptn_cod <- exons_prot_cod %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_name
  )) + labs(y = "Protein coding\ntranscripts") +
  geom_range(
    fill = "white",
    height = 0.25
  ) +
  geom_range(
    data = cds
  ) +
  geom_intron(
    data = to_intron(exons_prot_cod, "transcript_name"),
    aes(strand = strand),
    arrow.min.intron.length = 500 ) +
  ggtitle(title) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
    axis.line = element_line(colour = "black")) 



# junction data 
melted_intr_usage <- fread(paste0(outPath, "data4plots/adar_junctionTx_data.tsv"))
melted_intr_usage[ ,intron_usage:=round(intron_usage, 2)]

#plot for comparison
cols2Sel <- c( "chr", "start" , "end",  "strand", "intron_usage",  "transcript_name", "treatment")
comparisons <- unique(melted_intr_usage$comparison)

#post_fix <- "_ggtranscipt_colScheme_scaleSize.png"
#post_fix <- "_ggtranscipt_colScheme.png"


#######################
#Generate the plots
######################

library(patchwork)
#single plot - with space for the p110 and p150 ptn to be included
#####################################################################
ret <-lapply(comparisons, function(cmp){
  #cmp <- comparisons[1]
  l_dt <- melted_intr_usage[ comparison== cmp, ]
  setorder(l_dt, -treatment)
  tx <- sort(unique(l_dt$transcript_name))
  plots_tx_junct <-lapply(tx, function(txname){
    #txname <- tx[2]
    
    l_junct <- as_tibble(l_dt[ transcript_name == txname , .SD, .SDcols= cols2Sel ])
    col_pal <- palette_culture[ grep(paste0(l_junct$treatment, collapse= "|"), names(palette_culture))]
    col_pal <- col_pal[order(names(col_pal), decreasing = T)]
    # extract exons and cds for the MANE-select transcript
    l_tx_exons <- exons_prot_cod %>% dplyr::filter(transcript_name == txname)
    l_tx_cds <- cds %>% dplyr::filter(transcript_name == txname)
    return(generate_geom_junction_tx_plot(l_tx_exons, l_tx_cds, l_junct, col_pal ))
      
    
  })
  o_file <- paste0(outPath, unique(l_dt$culture), "_ggtranscript.png")
  
  plots <- plot_utr_ptn_cod +  plots_tx_junct[[1]] + plot_spacer() + plot_spacer()   + plots_tx_junct[[2]] + plot_spacer()  + plot_layout(nrow = 7, ncol = 1) +
    plot_annotation(tag_levels = 'i') 
  ggsave(o_file, plots, width = 8, height = 8)
   
  
})


#Single X axis - for the junctions
##########################################
ret <-lapply(comparisons, function(cmp){
  #cmp <- comparisons[1]
  l_dt <- melted_intr_usage[ comparison== cmp, ]
  
  setorder(l_dt, -transcript_name, -treatment)
  tx <- sort(unique(l_dt$transcript_name))
  l_junct <- as_tibble(l_dt[ , .SD, .SDcols= cols2Sel ])
  col_pal <- palette_culture[ grep(paste0(l_junct$treatment, collapse= "|"), names(palette_culture))]
  col_pal <- col_pal[order(names(col_pal), decreasing = T)]
  
  # extract exons and cds for the MANE-select transcript
  l_tx_exons <- exons_prot_cod %>% dplyr::filter(transcript_name %in% tx)
  l_tx_cds <- cds %>% dplyr::filter(transcript_name %in% tx)
  col_pal <- c(col_pal, col_pal)
  tib_exons<- l_tx_exons; tib_cds <- l_tx_cds; tib_junct <- l_junct
  rm( l_tx_exons, l_tx_cds, l_junct)
  p <- tib_exons %>%
             ggplot(aes(
               xstart = start,
               xend = end,
               y = transcript_name
             )) + labs(y = "") +
             geom_range(
               fill = "white", 
               height = 0.25
             ) +
             geom_range(
               data = tib_cds
             ) + 
             geom_intron(
               data = to_intron(tib_exons, "transcript_name"), strand="-"
             ) + 
             geom_junction(
               data = tib_junct,
               aes(size = intron_usage),
               junction.y.max = 0.5,
               show.legend = F, 
               color = col_pal 
             ) +
             geom_junction_label_repel(
               data = tib_junct,
               aes(label = round(intron_usage, 2)),
               junction.y.max = 0.5
             ) +
             theme(   axis.line = element_line(colour = "black") )

    
  o_file <- paste0(outPath, unique(l_dt$culture), "_ggtranscript_commonx.png")
  
  plots <- plot_utr_ptn_cod + plot_spacer() + p+ plot_layout(nrow = 3, ncol = 1, heights = c(1, 0.1, 1)) +
    plot_annotation(tag_levels = 'i') 
  ggsave(o_file, plots, width = 8, height = 9)
  
  
})

#excluding the protein coding trascript panel and including the 208 with the junctions
tx <- unique(exons_prot_cod$transcript_name)
ret <-lapply(comparisons, function(cmp){
  #cmp <- comparisons[1]
  l_dt <- melted_intr_usage[ comparison== cmp, ]
  
  setorder(l_dt, -transcript_name, -treatment)
  
  l_junct <- as_tibble(l_dt[ , .SD, .SDcols= cols2Sel ])
  col_pal <- palette_culture[ grep(paste0(l_junct$treatment, collapse= "|"), names(palette_culture))]
  col_pal <- col_pal[order(names(col_pal), decreasing = T)]
  
  # extract exons and cds for the MANE-select transcript
  l_tx_exons <- exons_prot_cod %>% dplyr::filter(transcript_name %in% tx)
  l_tx_cds <- cds %>% dplyr::filter(transcript_name %in% tx)
  col_pal <- c(col_pal, col_pal)
  tib_exons<- l_tx_exons; tib_cds <- l_tx_cds; tib_junct <- l_junct
  rm( l_tx_exons, l_tx_cds, l_junct)
  p <- tib_exons %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = transcript_name
    )) + labs(y = "") +
    geom_range(
      fill = "white", 
      height = 0.25
    ) +
    geom_range(
      data = tib_cds
    ) + 
    geom_intron(
      data = to_intron(tib_exons, "transcript_name"), strand="-"
    ) + 
    geom_junction(
      data = tib_junct,
      aes(size = intron_usage),
      junction.y.max = 0.5,
      show.legend = F, 
      color = col_pal 
    ) +
    geom_junction_label_repel(
      data = tib_junct,
      aes(label = round(intron_usage, 2)),
      junction.y.max = 0.5
    ) +
    theme(   axis.line = element_line(colour = "black") )
  
  
  o_file <- paste0(outPath, unique(l_dt$culture), "_ggtranscript_commonx_singlepanel.png")
  
  #ggsave(o_file, p, width = 8, height = 7)
  
  o_file <- paste0(outPath, unique(l_dt$culture), "_ggtranscript_commonx_singlepanel.pdf")
  
  ggsave(o_file, p, width = 8, height = 7)
  
})



#generating individual transcript plots for the astrocyte and coculture 
##########################################################################
img_width <- 9
img_height <- 4
ggsave( paste0(outPath, "ggtranscript_utr_ptn_cod.png"),plot_utr_ptn_cod, width = img_width, height = img_height, dpi=300, units = "in")

#plot for comparison
cols2Sel <- c( "chr", "start" , "end",  "strand", "intron_usage",  "transcript_name", "treatment")
comparisons <- unique(melted_intr_usage$comparison)

ret <-lapply(comparisons, function(cmp){
  #cmp <- comparisons[1]
  l_dt <- melted_intr_usage[ comparison== cmp, ]
  setorder(l_dt, -treatment)
  tx <- sort(unique(l_dt$transcript_name))
  plots_tx_junct <-lapply(tx, function(txname){
        l_junct <- as_tibble(l_dt[ transcript_name == txname , .SD, .SDcols= cols2Sel ])
        col_pal <- palette_culture[ grep(paste0(l_junct$treatment, collapse= "|"), names(palette_culture))]
        col_pal <- col_pal[order(names(col_pal), decreasing = T)]
        # extract exons and cds for the MANE-select transcript
        l_tx_exons <- exons_prot_cod %>% dplyr::filter(transcript_name == txname)
        l_tx_cds <- cds %>% dplyr::filter(transcript_name == txname)
        return(generate_geom_junction_tx_plot(l_tx_exons, l_tx_cds, l_junct, col_pal ))

    })

  ggsave( paste0(outPath, "ggtranscript_juncTx_", tx[1] ,"_",unique(l_dt$culture),".png"),plots_tx_junct[[1]], width = img_width, height = img_height, dpi=300, units = "in")
  ggsave( paste0(outPath, "ggtranscript_juncTx_", tx[2] ,"_",unique(l_dt$culture),".png"),plots_tx_junct[[2]], width = img_width, height = img_height, dpi=300, units = "in")
})


generate_geom_junction_tx_plot <- function(tib_exons, tib_cds, tib_junct, col_pal){
  #tib_exons<- l_tx_exons; tib_cds <- l_tx_cds; tib_junct <- l_junct
  return(tib_exons %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = transcript_name
    )) + labs(y = "") +
    geom_range(
      fill = "white", 
      height = 0.25
    ) +
    geom_range(
      data = tib_cds
    ) + 
    geom_intron(
      data = to_intron(tib_exons, "transcript_name"), strand="-"
    ) + 
    geom_junction(
      data = tib_junct,
      aes(size = intron_usage),
      junction.y.max = 0.5,
      show.legend = F, 
      color = col_pal 
    ) +
    geom_junction_label_repel(
      data = tib_junct,
      aes(label = round(intron_usage, 2)),
      junction.y.max = 0.5
    ) +
      theme(   axis.line = element_line(colour = "black"))# + 
    #scale_size_continuous(range = c(1, 5)) #c(0.1, 1)
)
}


