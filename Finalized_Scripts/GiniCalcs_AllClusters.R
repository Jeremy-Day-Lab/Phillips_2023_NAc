#Calculate gini coefficients for all clusters
library(reldist)
library(Seurat)
library(ggplot2)
library(scales)
library(ggrepel)

#Load the NAc_Combo object
NAc_Combo <- readRDS(file = "2019-JD-0040/MCN_Code/Objects/NAc_Combo_Integrated.RDS")

DimPlot(object = NAc_Combo,reduction = "umap",label = TRUE) + NoLegend()

#Set the Seurat assay to RNA
DefaultAssay(NAc_Combo) <- "RNA"

#Set identities back to the integer values instead of the characters used for DEG testing
Idents(NAc_Combo) <- NAc_Combo$integrated_snn_res.0.15

#Calculate the cluster-specific average expression values for every gene
CellAverages <- AverageExpression(NAc_Combo)

#Pull the RNA values
RNA <- CellAverages$RNA

#Calculate the Gini coefficient, maximum/mean gene expression value, and ratio
Gini_calc <- apply(RNA, 1, gini)
max_RNA <- apply(RNA, 1, max)
mean_RNA <- apply(RNA, 1, mean)
Ratio <- max_RNA/mean_RNA

#Create a dataframe consisting of the Gini coefficient, maximum/mean gene expression value, and ratio
Gini_dataframe <- data.frame(Gini_calc, max_RNA, mean_RNA, Ratio, RNA)
#Create a column of gene names
Gini_dataframe$Gene <- rownames(Gini_dataframe)

# "0" = "Olig-1",#
# "1" = "Drd2-MSN-1",#
# "2" = "Drd1-MSN-1",#
# "3" = "GABAergic",#
# "4" = "Grm8-MSN",#
# "5" = "Astrocyte",#
# "6" = "Drd1-MSN-2",#
# "7" = "Polydendrocyte",
# "8" = "Microglia",#
# "9" = "Pvalb-Interneuron",#
# "10" = "Glutamatergic",#
# "11" = "Drd3-MSN",#
# "12" = "Drd2-MSN-2",#
# "13" = "Sst-Interneuron",#
# "14" = "Mural",#
# "15" = "Chat-Interneuron"
colnames(Gini_dataframe)[5:20] <- c("Olig-1","Drd2-MSN-1","Drd1-MSN-1","GABAergic",
                                    "Grm8-MSN","Astrocyte","Drd1-MSN-2","Polydendrocyte",
                                    "Microglia","Pvalb-Interneuron","Glutamatergic","Drd3-MSN",
                                    "Drd2-MSN-2","Sst-Interneuron","Mural","Chat-Interneuron")

#Remove any NA Gini_calc values
Gini_dataframe_noNA <- Gini_dataframe[-which(is.na(Gini_dataframe$Gini_calc)),]

#Reorder the column names
Gini_dataframe_noNA <- Gini_dataframe_noNA[,c("Gini_calc","max_RNA","mean_RNA","Ratio",
                                              "Polydendrocyte","Olig-1","Mural","Microglia",
                                              "Astrocyte","Glutamatergic","Sst-Interneuron","Pvalb-Interneuron",
                                              "Chat-Interneuron","GABAergic","Grm8-MSN","Drd3-MSN",
                                              "Drd2-MSN-2","Drd2-MSN-1","Drd1-MSN-2","Drd1-MSN-1","Gene")]


#Change all of the - to .
colnames(Gini_dataframe_noNA)[grep(pattern = "-",x = colnames(Gini_dataframe_noNA))] <- gsub(colnames(Gini_dataframe_noNA)[grep(pattern = "-",x = colnames(Gini_dataframe_noNA))],pattern = "-",replacement = ".")

#Write out the gini data.frame
write.table(x         = Gini_dataframe_noNA,
            file      = "2019-JD-0040/MCN_Code/Tables/Gini_Dataframe_AllClusters.txt",
            col.names = TRUE,
            row.names = TRUE,
            sep       = "\t",
            quote     = FALSE)
 

Idents(NAc_Combo) <- NAc_Combo$Combo_CellType
#Make Gini plots
for(i in colnames(Gini_dataframe_noNA)[5:20]){
  print(i)
  #Pull the DEG list
  x <- read.table(file = paste0("2019-JD-0040/MCN_Code/Tables/DESeq2_CellTypes/",i,".txt"), sep = "\t",header = TRUE)
  #Calculate the percent of cells expressing 
  celltype_name <- gsub(pattern = "[.]",replacement = "-",x = i)
  cells.1    <- WhichCells(NAc_Combo,idents = celltype_name)
  counts_mat <- NAc_Combo@assays$RNA@counts
  x$Pct_Expressing <- NA
  x$Pct_Expressing <- (rowSums(x = counts_mat[x$GeneName, cells.1,drop = FALSE] > 0) / length(cells.1))*100
  #Merge the dataframes 
  y <- merge(x = x,
             y = Gini_dataframe_noNA,
             by.x = "GeneName",
             by.y = "Gene")
  y$logCluster <- log10(y[,i])
  y <- y[!is.infinite(y$logCluster),]
  #Make the oplot
  Gini_plot <- ggplot(data = y,aes(x = logCluster,y = Gini_calc,size = Pct_Expressing)) +
    geom_point(color = "lightgrey",aes(size = Pct_Expressing)) +
    geom_point(data = subset(y,subset=(padj <= 0.05 & log2FoldChange > 0)),
               color = hue_pal()(16)[which(colnames(Gini_dataframe_noNA)[5:20] == i)]) +
    geom_text_repel(data = subset(y,subset=(padj <= 0.05 & log2FoldChange > 0))[order(subset(y,subset=(padj <= 0.05 & log2FoldChange > 0))$Gini_calc,decreasing = TRUE)[1:10],],
                    label = subset(y,subset=(padj <= 0.05 & log2FoldChange > 0))[order(subset(y,subset=(padj <= 0.05 & log2FoldChange > 0))$Gini_calc,decreasing = TRUE)[1:10],"GeneName"],
                    show.legend = FALSE,
                    size = 6) +
    geom_density2d(color = "red",show.legend = FALSE) +
    ylim(c(0,1)) +
    theme_bw() + 
    NoGrid() +
    labs(x = paste0("log10(",i," Counts)"),
         y = "Gini coefficient",
         size = "% of Cells in Cluster\nExpressing Gene") +
    ggtitle(label = i) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(plot     = Gini_plot,
         filename = paste0("/data/project/daylab/2019-JD-0040/MCN_Code/Plots/DESeq2_DEGResults/Gini_Plots/",i,".pdf"),
         height   = 12,
         width    = 12)
  print(paste(i,"is complete"))
}



sessionInfo()
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 16.04.6 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/openblas-base/libblas.so.3
# LAPACK: /usr/lib/libopenblasp-r0.2.18.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8      
# [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggrepel_0.8.2      scales_1.1.1       ggplot2_3.3.2      SeuratObject_4.0.2 Seurat_4.0.4       reldist_1.6-6     
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15            colorspace_1.4-1      deldir_1.0-6          ellipsis_0.3.2        ggridges_0.5.2        htmlTable_2.1.0       base64enc_0.1-3       rstudioapi_0.11       spatstat.data_2.1-0  
# [10] farver_2.0.3          leiden_0.3.3          listenv_0.8.0         codetools_0.2-16      splines_4.0.2         knitr_1.30            polyclip_1.10-0       Formula_1.2-4         jsonlite_1.7.1       
# [19] ica_1.0-2             cluster_2.1.0         png_0.1-7             uwot_0.1.10           shiny_1.5.0           sctransform_0.3.2     spatstat.sparse_2.0-0 compiler_4.0.2        httr_1.4.2           
# [28] backports_1.1.10      Matrix_1.3-4          fastmap_1.1.0         lazyeval_0.2.2        cli_3.3.0             later_1.1.0.1         htmltools_0.5.2       tools_4.0.2           igraph_1.2.6         
# [37] gtable_0.3.0          glue_1.6.2            RANN_2.6.1            reshape2_1.4.4        dplyr_1.0.2           rappdirs_0.3.1        Rcpp_1.0.7            scattermore_0.7       vctrs_0.4.1          
# [46] nlme_3.1-149          lmtest_0.9-38         xfun_0.18             stringr_1.4.0         globals_0.13.1        mime_0.9              miniUI_0.1.1.1        lifecycle_0.2.0       irlba_2.3.3          
# [55] goftest_1.2-2         future_1.19.1         MASS_7.3-53           zoo_1.8-8             spatstat.core_2.3-0   promises_1.1.1        spatstat.utils_2.2-0  parallel_4.0.2        RColorBrewer_1.1-2   
# [64] yaml_2.2.1            reticulate_1.16       pbapply_1.4-3         gridExtra_2.3         rpart_4.1-15          latticeExtra_0.6-29   stringi_1.7.6         checkmate_2.0.0       rlang_1.0.2          
# [73] pkgconfig_2.0.3       matrixStats_0.57.0    lattice_0.20-41       ROCR_1.0-11           purrr_0.3.4           tensor_1.5            labeling_0.3          patchwork_1.0.1       htmlwidgets_1.5.2    
# [82] cowplot_1.1.0         tidyselect_1.1.0      RcppAnnoy_0.0.19      plyr_1.8.6            magrittr_1.5          R6_2.4.1              generics_0.0.2        Hmisc_4.4-1           withr_2.3.0          
# [91] pillar_1.4.6          foreign_0.8-80        mgcv_1.8-33           fitdistrplus_1.1-1    survival_3.2-7        abind_1.4-5           nnet_7.3-14           tibble_3.0.4          future.apply_1.6.0   
# [100] crayon_1.3.4          KernSmooth_2.23-17    spatstat.geom_2.4-0   plotly_4.9.2.1        jpeg_0.1-8.1          isoband_0.2.2         grid_4.0.2            data.table_1.13.0     digest_0.6.26        
# [109] xtable_1.8-4          tidyr_1.1.2           httpuv_1.5.4          munsell_0.5.0         viridisLite_0.3.0    
