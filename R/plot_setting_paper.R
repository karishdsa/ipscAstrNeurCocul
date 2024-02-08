library(ggsci)
library(scales)

#variables
SIZE_LINE_YINT<-0.25
SIZE_GEOM_LINE<-0.5
SIZE_GEOM_POINT <- 1 #0.75
SIZE_AXIS_TITLE <-8 #7
SIZE_AXIS_TEXT <-8 #7
SIZE_LEGEND_TITLE <-8 #7
SIZE_LEGEND_TEXT <- 8#7
SIZE_PLOT_TITLE <- 9

JITTER_WIDTH <- 0.2

IMAGE_DPI_COLOR <- 300
IMAGE_DPI_BLACKANDWHITE <- 500

IMAGE_WIDTH_1COLUMN <- 85
IMAGE_WIDTH_1.5COLUMN <- 114
IMAGE_WIDTH_2COLUMN <- 174
IMAGE_WIDTH_UNITS <- "mm"

BOXPLOT_WIDTH <- 0.5

#colours
COLOUR_CUTOFF_RED <- pal_aaas()(10)[2] #red
ALPHA_VAL <- 0.7


palette_clusters <- c( "Astrocyte cluster 1_untreated" =  "#B2FFFF", 
                       "Astrocyte cluster 2_untreated" = "#55E4F2",
                       "Neuron cluster 1_untreated" = "#ffe4ff", 
                       "Neuron cluster 2_untreated" = "#ffa3ff", 
                       "Neuron cluster 3_untreated" = "#EA60C2",  
                       "Astrocyte cluster 1_treated" =  "#85BFBF",
                       "Astrocyte cluster 2_treated" = "#3b9fa9",
                       "Neuron cluster 1_treated" = "#e5a0e5", 
                       "Neuron cluster 2_treated" = "#cc51cc", 
                       "Neuron cluster 3_treated" = "#a34387"
                       
)
palette_culture <-c("neuron_treated" =  "#CD1076",
                    "neuron_untreated"=   "#FF91C8",
                    "astrocyte_treated"= "#0097A6",
                    "astrocyte_untreated"= "#7FDEEA",
                    "coculture_treated"= "#303F9F",
                    "coculture_untreated"= "#9FA7D8" )

palette_clusters_comb <- c( "Astrocyte cluster 1" =  "#B2FFFF", 
                            "Astrocyte cluster 2" = "#55E4F2", 
                            "Neuron cluster 1" = "#FFB2FF", 
                            "Neuron cluster 2" = "#FF66FF",
                            "Neuron cluster 3" = "#EA60C2", 
                            "Neuron cluster 4" = "#FF0065", 
                            "Neuron cluster 5" = "#FB7979", 
                            "Neuron cluster 6" = "#B000B0")

#set the theme with base font
theme_plot <- function( base_family = "Arial") {
  theme_bw( base_family = base_family ) +
    theme(legend.position = "bottom") +
    theme(panel.grid.major = element_blank(), 
          #      panel.background = element_blank(), 
          panel.border = element_blank(),
          #       axis.line = element_line(colour = "black"), #clear the background
          strip.background = element_blank()) +
    theme(axis.title  = element_text(face="italic", 
                                     color="black", 
                                     size=SIZE_AXIS_TITLE),
          axis.text  = element_text( color="black", 
                                     size=SIZE_AXIS_TEXT) ) +
    theme(plot.title = element_text(size=SIZE_PLOT_TITLE, 
                                    hjust=0.5)) +
    theme( legend.title = element_text(colour="black", 
                                       size=SIZE_LEGEND_TITLE), 
           legend.text = element_text(colour="black", 
                                      size=SIZE_LEGEND_TEXT) ) +
    theme( plot.tag = element_text(size = SIZE_LEGEND_TITLE) )+
    theme(plot.subtitle=element_text(size=SIZE_PLOT_TITLE, 
                                     hjust=0.5))
}

#update the font defaults in the plots
update_font_defaults <- function(font ="Arial") {
  
  update_geom_defaults("text", list(family = font))
  update_geom_defaults("label", list(family = font))
  
}

#set the theme 
set_plot_theme <- function(font="Arial")  {
  
  theme_set(theme_plot(font))
  
  update_font_defaults( font)
  
}

set_plot_theme("Helvetica")

