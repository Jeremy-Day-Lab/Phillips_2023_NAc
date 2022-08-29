#Goal for this analysis is to integrate the repeated and acute objects

#Set working directory and seed
setwd("/data/project/daylab/")
set.seed(1234)

#Load seurat and ggplot2
library(Seurat)
library(ggplot2)

#Load the repeated object
Repeated <- readRDS(file = "2019-JD-0040/MCN_Code/Objects/Repeated_Cocaine_DubsIdentified.RDS")

#Make the dataset metadata 
Repeated$Dataset <- "Repeated"

#Now subset for singlets only
#26220 cells to start
Repeated <- subset(Repeated,subset=(Doublet_Classification == "Singlet"))
#23599 cells post subset

#Now load the acute object
Acute <- readRDS(file = "/data/project/daylab/2019-JD-0040/MCN_Code/Objects/Acute_Decontaminated.RDS")
#15655 cells
#Acute had no doublets to subset or anything. 

#Now integrate both datasets. 
#Move everything together 
NAc_Combo <- FindIntegrationAnchors(object.list = list(Acute,Repeated), dims = 1:17)
NAc_Combo  <- IntegrateData(anchorset = NAc_Combo,dims = 1:17)

DefaultAssay(NAc_Combo) <- "integrated"

# Run the standard workflow for visualization and clustering
NAc_Combo <- ScaleData(NAc_Combo,verbose = FALSE)
NAc_Combo <- RunPCA(NAc_Combo,npcs = 17 ,verbose = FALSE) #Compute 50 npcs by default
# Dimensionality reduction and Clustering
NAc_Combo <- RunUMAP(NAc_Combo, reduction = "pca", dims = 1:17)
NAc_Combo <- FindNeighbors(NAc_Combo, reduction = "pca", dims = 1:17)

NAc_Combo <- FindClusters(NAc_Combo, resolution = 0.15)

DimPlot(NAc_Combo,reduction = "umap",label = TRUE) + NoLegend()

saveRDS(object = NAc_Combo,file = "2019-JD-0040/MCN_Code/Objects/NAc_Combo_v0.RDS")

Dp1 <- DotPlot(object = NAc_Combo,
               assay  = "RNA",
               features =  c("Drd1","Pdyn","Ebf1", #D1-MSN
                             "Drd2","Penk", #D2-MSN
                             "Drd3","Grm8", #D3/Grm8 MSN
                             "Elavl2","Kit","Sst","Chat",#GABAergic markers
                             "Slc17a7", #Glut
                             "Mbp","Opalin", #Oligs
                             "Aqp4","Gja1", #Astrocytes
                             "Pdgfra", #Polydendrocytes
                             "Arhgap15", #Microglia
                             "Rgs5", #Mural
                             "Ppp1r1b","Foxp2","Bcl11b","Gad1","Syt1"), #Neuronal markers
               cols = c("lightgrey","red")) +
  theme(axis.text.x = element_text(angle = 45, hjust =1))


#0: Olig-1
#1: Drd2-MSN-1
#2: Drd1-MSN-1
#3: GABAergic-1
#4: Grm8-MSN
#5: Astrocyte
#6: Drd1-MSN-2
#7: Polydendrcyte
#8: Microglia
#9: Pvalb-Interneuron
#10: Glutamatergic
#11: Drd3-MSN
#12: Drd2-MSN-2
#13: Sst-Interneuron
#14: Mural
#15: Chat-Interneuron

#Change the identities
NAc_Combo <-  RenameIdents(object = NAc_Combo,
                           "0" = "Olig-1",#
                           "1" = "Drd2-MSN-1",#
                           "2" = "Drd1-MSN-1",#
                           "3" = "GABAergic",#
                           "4" = "Grm8-MSN",#
                           "5" = "Astrocyte",#
                           "6" = "Drd1-MSN-2",#
                           "7" = "Polydendrocyte",
                           "8" = "Microglia",#
                           "9" = "Pvalb-Interneuron",#
                           "10" = "Glutamatergic",#
                           "11" = "Drd3-MSN",#
                           "12" = "Drd2-MSN-2",#
                           "13" = "Sst-Interneuron",#
                           "14" = "Mural",#
                           "15" = "Chat-Interneuron")#


#Change the order of the identities
Idents(NAc_Combo) <- factor(x = Idents(NAc_Combo),
                            levels = rev(c("Drd1-MSN-1",
                                           "Drd1-MSN-2",
                                           "Drd2-MSN-1",
                                           "Drd2-MSN-2",
                                           "Drd3-MSN",
                                           "Grm8-MSN",
                                           "GABAergic",
                                           "Chat-Interneuron",
                                           "Pvalb-Interneuron",
                                           "Sst-Interneuron",
                                           "Glutamatergic",
                                           "Astrocyte",
                                           "Microglia",
                                           "Mural",
                                           "Olig-1",
                                           "Polydendrocyte")))

NAc_Combo$Combo_CellType <- Idents(NAc_Combo)


#remake the dotplot
DefaultAssay(NAc_Combo) <- "RNA"
Dp1 <- DotPlot(object = NAc_Combo,
               features =  c("Drd1","Pdyn","Ebf1", #D1-MSN
                             "Drd2","Penk", #D2-MSN
                             "Drd3","Grm8", #D3/Grm8 MSN
                             "Elavl2","Chat","Kit","Sst",#GABAergic markers
                             "Slc17a7", #Glut
                             "Mbp","Opalin", #Oligs
                             "Aqp4","Gja1", #Astrocytes
                             "Pdgfra", #Polydendrocytes
                             "Arhgap15", #Microglia
                             "Rgs5", #Mural
                             "Ppp1r1b","Foxp2","Bcl11b","Gad1","Syt1"), #Neuronal markers
               cols = c("lightgrey","red"),
               dot.scale = 13) +
  theme(axis.text.x = element_text(angle = 90,vjust = 0))

#Save the dotplot
ggsave(plot = Dp1,
       filename = "2019-JD-0040/MCN_Code/Plots/Integrated_Dotplot_with_CellTypeLabels.pdf",
       height = 12,
       width = 12)

Idents(NAc_Combo) <- factor(x = Idents(NAc_Combo),
                            levels = rev(levels(Idents(NAc_Combo))))
VP1 <- VlnPlot(object = NAc_Combo,
        features =  c("Drd1","Ebf1", #D1-MSN
                      "Drd2","Penk", #D2-MSN
                      "Drd3","Grm8", #D3/Grm8 MSN
                      "Elavl2","Chat","Kit","Sst",#GABAergic markers
                      "Slc17a7", #Glut
                      "Aqp4","Gja1", #Astrocytes
                      "Arhgap15", #Microglia
                      "Rgs5", #Mural
                      "Mbp","Opalin", #Oligs
                      "Pdgfra", #Polydendrocytes,
                      "Ppp1r1b","Foxp2","Bcl11b","Gad1","Syt1"),
        flip = TRUE,
        stack = TRUE,
        fill.by = "ident") + 
  NoLegend()
ggsave(plot = VP1,
       filename = "2019-JD-0040/MCN_Code/Plots/Integrated_VlnPlot.pdf",
       height = 7,
       width = 3.5)

#Save some featureplots
#Syt1
Syt1 <- FeaturePlot(object = NAc_Combo,features = "Syt1",cols = c("lightgrey","red"),max.cutoff = 4) + NoAxes()
ggsave(plot = Syt1,
       filename = "2019-JD-0040/MCN_Code/Plots/Syt1_FeaturePlot_Integrated.pdf",
       height = 8,
       width  = 8)
ggsave(plot = Syt1 + ggtitle("") + NoLegend(),
       filename = "2019-JD-0040/MCN_Code/Plots/Syt1_FeaturePlot_Integrated_NoTitle.jpg",
       height = 8,
       width  = 8)
#Gad1
Gad1 <- FeaturePlot(object = NAc_Combo,features = "Gad1",cols = c("lightgrey","red"),max.cutoff = 4) + NoAxes()
ggsave(plot = Gad1,
       filename = "2019-JD-0040/MCN_Code/Plots/Gad1_FeaturePlot_Integrated.pdf",
       height = 8,
       width  = 8)
ggsave(plot = Gad1 + ggtitle("") + NoLegend(),
       filename = "2019-JD-0040/MCN_Code/Plots/Gad1_FeaturePlot_Integrated_NoTitle.jpg",
       height = 8,
       width  = 8)
#Bcl11b 
Bcl11b <- FeaturePlot(object = NAc_Combo,features = "Bcl11b",cols = c("lightgrey","red"),max.cutoff = 4) + NoAxes()
ggsave(plot = Bcl11b,
       filename = "2019-JD-0040/MCN_Code/Plots/Bcl11b_FeaturePlot_Integrated.pdf",
       height = 8,
       width  = 8)
ggsave(plot = Bcl11b + ggtitle("") + NoLegend(),
       filename = "2019-JD-0040/MCN_Code/Plots/Bcl11b_FeaturePlot_Integrated_NoTitle.jpg",
       height = 8,
       width  = 8)
#Arhgap15 
Arhgap15 <- FeaturePlot(object = NAc_Combo,features = "Arhgap15",cols = c("lightgrey","red"),max.cutoff = 4) + NoAxes()
ggsave(plot = Arhgap15,
       filename = "2019-JD-0040/MCN_Code/Plots/Arhgap15_FeaturePlot_Integrated.pdf",
       height = 8,
       width  = 8)
ggsave(plot = Arhgap15 + ggtitle("") + NoLegend(),
       filename = "2019-JD-0040/MCN_Code/Plots/Arhgap15_FeaturePlot_Integrated_NoTitle.jpg",
       height = 8,
       width  = 8)
#Mbp 
Mbp <- FeaturePlot(object = NAc_Combo,features = "Mbp",cols = c("lightgrey","red"),max.cutoff = 4) + NoAxes()
ggsave(plot = Mbp,
       filename = "2019-JD-0040/MCN_Code/Plots/Mbp_FeaturePlot_Integrated.pdf",
       height = 8,
       width  = 8)
ggsave(plot = Mbp + ggtitle("") + NoLegend(),
       filename = "2019-JD-0040/MCN_Code/Plots/Mbp_FeaturePlot_Integrated_NoTitle.jpg",
       height = 8,
       width  = 8)
#Gja1 
Gja1 <- FeaturePlot(object = NAc_Combo,features = "Gja1",cols = c("lightgrey","red"),max.cutoff = 4) + NoAxes()
ggsave(plot = Gja1,
       filename = "2019-JD-0040/MCN_Code/Plots/Gja1_FeaturePlot_Integrated.pdf",
       height = 8,
       width  = 8)
ggsave(plot = Gja1 + ggtitle("") + NoLegend(),
       filename = "2019-JD-0040/MCN_Code/Plots/Gja1_FeaturePlot_Integrated_NoTitle.jpg",
       height = 8,
       width  = 8)

#Make some umaps
CellType_umap <- DimPlot(object = NAc_Combo,label = TRUE) + NoLegend() + NoAxes() + ggtitle("Joint") + theme(plot.title = element_text(hjust = 0.5))
Acute_umap <- DimPlot(object = subset(NAc_Combo,subset=(Dataset == "Acute")),label = TRUE) + NoLegend() + NoAxes() + ggtitle("Acute") + theme(plot.title = element_text(hjust = 0.5))
Repeated_umap <- DimPlot(object = subset(NAc_Combo,subset=(Dataset == "Repeated")),label = TRUE) + NoLegend() + NoAxes() + ggtitle("Repeated") + theme(plot.title = element_text(hjust = 0.5))
cowplot::plot_grid(plotlist = list(CellType_umap,Acute_umap,Repeated_umap),ncol = 3)
#Save these plots
ggsave(plot = CellType_umap,
       filename = "2019-JD-0040/MCN_Code/Plots/Integrated_Umap_with_CellTypeLabels_NoAxes.jpg",
       height = 12,
       width = 12)
ggsave(plot = Acute_umap,
       filename = "2019-JD-0040/MCN_Code/Plots/Acute_Umap_with_CellTypeLabels_NoAxes.jpg",
       height = 12,
       width = 12)
ggsave(plot = Repeated_umap,
       filename = "2019-JD-0040/MCN_Code/Plots/Repeated_Umap_with_CellTypeLabels_NoAxes.jpg",
       height = 12,
       width = 12)
#Make umaps without labels
CellType_umap <- DimPlot(object = NAc_Combo,label = FALSE) + NoLegend() + NoAxes() + ggtitle("Joint") + theme(plot.title = element_text(hjust = 0.5))
Acute_umap <- DimPlot(object = subset(NAc_Combo,subset=(Dataset == "Acute")),label = FALSE) + NoLegend() + NoAxes() + ggtitle("Acute") + theme(plot.title = element_text(hjust = 0.5))
Repeated_umap <- DimPlot(object = subset(NAc_Combo,subset=(Dataset == "Repeated")),label = FALSE) + NoLegend() + NoAxes() + ggtitle("Repeated") + theme(plot.title = element_text(hjust = 0.5))
#Save these plots
ggsave(plot = CellType_umap,
       filename = "2019-JD-0040/MCN_Code/Plots/Integrated_Umap_without_CellTypeLabels_NoAxes.jpg",
       height = 12,
       width = 12)
ggsave(plot = Acute_umap,
       filename = "2019-JD-0040/MCN_Code/Plots/Acute_Umap_without_CellTypeLabels_NoAxes.jpg",
       height = 12,
       width = 12)
ggsave(plot = Repeated_umap,
       filename = "2019-JD-0040/MCN_Code/Plots/Repeated_Umap_without_CellTypeLabels_NoAxes.jpg",
       height = 12,
       width = 12)

#Make a umap colored by dataset
Idents(NAc_Combo) <- NAc_Combo$Dataset
Umap_Labled_by_dataset <- DimPlot(object = NAc_Combo,label = FALSE) + NoAxes() + ggtitle("Joint") + theme(plot.title = element_text(hjust = 0.5))
ggsave(plot = Umap_Labled_by_dataset,
       filename = "2019-JD-0040/MCN_Code/Plots/Umap_Labled_by_dataset_NoAxes.jpg",
       height = 12,
       width = 12)

#make a table of the number of cells in each cluster
Idents(NAc_Combo) <- NAc_Combo$Combo_CellType
CellType_Freqs <- as.data.frame(table(Idents(NAc_Combo)))
row.names(CellType_Freqs) <- CellType_Freqs$Var1
CellType_Freqs$Prop_Repeated <- NA
for(i in CellType_Freqs$Var1){
  print(i)
  CellType_Freqs[i,"Prop_Repeated"] <- ncol(subset(NAc_Combo,subset=(Combo_CellType == i & Dataset == "Repeated")))/23599
  CellType_Freqs[i,"Prop_Acute"]    <- ncol(subset(NAc_Combo,subset=(Combo_CellType == i & Dataset == "Acute")))/15655
}
colnames(CellType_Freqs)[1:2] <- c("CellType","Frequency")
CellType_Freqs$Prop_Repeated <- CellType_Freqs$Prop_Repeated*100
CellType_Freqs$Prop_Acute    <- CellType_Freqs$Prop_Acute*100
colnames(CellType_Freqs)[3:4] <- c("Repeated","Acute")

CellType_Freqs_melt <- reshape2::melt(CellType_Freqs[,c("CellType","Repeated","Acute")])

bp1 <- ggplot(data = CellType_Freqs_melt,aes(y = reorder(CellType,value), x = value, fill = variable)) +
  geom_bar(stat = "identity",position = "dodge") +
  labs(y    = "Cell Type",
       x    = "% of Cells Captured",
       fill = "Dataset") +
  theme_bw() +
  NoGrid()

ggsave(plot = bp1, filename = "2019-JD-0040/MCN_Code/Plots/Barplot_prop_of_cells_captured.pdf",height = 12, width = 12)


###Save the integrated object
saveRDS(object = NAc_Combo,
        file = "2019-JD-0040/MCN_Code/Objects/NAc_Combo_Integrated.RDS")

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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_3.3.2      SeuratObject_4.0.2 Seurat_4.0.4      
# 
# loaded via a namespace (and not attached):
#   [1] nlme_3.1-149          matrixStats_0.57.0    spatstat.sparse_2.0-0 RcppAnnoy_0.0.19      RColorBrewer_1.1-2    httr_1.4.2            sctransform_0.3.2     tools_4.0.2          
# [9] R6_2.4.1              irlba_2.3.3           rpart_4.1-15          KernSmooth_2.23-17    uwot_0.1.10           mgcv_1.8-33           lazyeval_0.2.2        colorspace_1.4-1     
# [17] withr_2.3.0           tidyselect_1.1.0      gridExtra_2.3         compiler_4.0.2        cli_3.3.0             plotly_4.9.2.1        labeling_0.3          scales_1.1.1         
# [25] lmtest_0.9-38         spatstat.data_2.1-0   ggridges_0.5.2        pbapply_1.4-3         rappdirs_0.3.1        goftest_1.2-2         stringr_1.4.0         digest_0.6.26        
# [33] spatstat.utils_2.2-0  pkgconfig_2.0.3       htmltools_0.5.2       fastmap_1.1.0         htmlwidgets_1.5.2     rlang_1.0.2           rstudioapi_0.11       shiny_1.5.0          
# [41] farver_2.0.3          generics_0.0.2        zoo_1.8-8             jsonlite_1.7.1        ica_1.0-2             dplyr_1.0.2           magrittr_1.5          patchwork_1.0.1      
# [49] Matrix_1.3-4          Rcpp_1.0.7            munsell_0.5.0         abind_1.4-5           reticulate_1.16       lifecycle_0.2.0       stringi_1.7.6         yaml_2.2.1           
# [57] MASS_7.3-53           Rtsne_0.15            plyr_1.8.6            grid_4.0.2            parallel_4.0.2        listenv_0.8.0         promises_1.1.1        ggrepel_0.8.2        
# [65] crayon_1.3.4          miniUI_0.1.1.1        deldir_1.0-6          lattice_0.20-41       cowplot_1.1.0         splines_4.0.2         tensor_1.5            pillar_1.4.6         
# [73] igraph_1.2.6          spatstat.geom_2.4-0   future.apply_1.6.0    reshape2_1.4.4        codetools_0.2-16      leiden_0.3.3          glue_1.6.2            data.table_1.13.0    
# [81] png_0.1-7             vctrs_0.4.1           httpuv_1.5.4          gtable_0.3.0          RANN_2.6.1            purrr_0.3.4           spatstat.core_2.3-0   polyclip_1.10-0      
# [89] tidyr_1.1.2           scattermore_0.7       future_1.19.1         mime_0.9              xtable_1.8-4          RSpectra_0.16-0       later_1.1.0.1         survival_3.2-7       
# [97] viridisLite_0.3.0     tibble_3.0.4          cluster_2.1.0         globals_0.13.1        fitdistrplus_1.1-1    ellipsis_0.3.2        ROCR_1.0-11
