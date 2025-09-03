library('pheatmap')
library('RColorBrewer')
library(ggplot2)
var=read.csv("C:/Dropbox/My PC (Andrei-PC)/Desktop/mapping heatmap/RPMM/final matrix/all.txt", header = TRUE,sep = "\t")
#var=read.csv("C:/Dropbox/My PC (Andrei-PC)/Desktop/mapping heatmap/correlation/global/sorted_matrix1.txt", header = TRUE,sep = "\t")

rownames(var)=var[,1]
var[,1]=NULL
#test=test[-65,]

#all the value less than 1 becomes 1 for log2 treatment
var[var<1]=1
var=as.matrix(var)
var2=log2(var)

#check the distribution of data.


#define quantile breaks
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

#set how many breaks you want to give. You have to provide n-1 number of colors in the heatmap
mat_breaks <- quantile_breaks(var2, n = 12)

pheatmap(
  mat               = var2,
  color             = rev(brewer.pal(11,"Spectral")),
  breaks            = mat_breaks,
  border_color      = NA,
  show_colnames     = T,
  show_rownames     = T,
  drop_levels       = TRUE,
  fontsize          = 14,
  cluster_rows = T, cluster_cols = F
  ,
  angle_col         = 45 
)

#mat_breaks <- quantile_breaks(test2, n = 30)
#install.packages("paletteer")

#library(paletteer)
#library(ggthemes)
#color_palette <- paletteer_c("grDevices::Plasma", 30)


#pheatmap(
  #mat               = test2,
  #color             = rev(brewer.pal(11, "Spectral")),
  #color             = color_palette,
  #breaks            = mat_breaks,
  #border_color      = NA,
  #show_colnames     = TRUE,
  #show_rownames     = TRUE,
  #drop_levels       = TRUE,
  #fontsize          = 14,
  #cluster_rows      = TRUE,
  #cluster_cols      = FALSE,
  #angle_col         = 45  # Rotate column headers by 45 degrees
#)
