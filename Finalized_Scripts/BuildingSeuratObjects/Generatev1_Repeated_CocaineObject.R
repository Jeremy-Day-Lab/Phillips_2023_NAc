#Data from this analysis comes from male and female adult Sprague-Dawley rats that were treated for 7 consecutive days 
#with 20mg/kg of cocaine. 1 hour forlloiwng final exposure, animals were sacrificed, and NAc was taken. 
#Goal of this analysis: Generate a seurat object for this dataset. We will later use SoupX to remove transcripts found overrepresented
#in empty droplets, and doubletfinder to remove doublets. 

#Set working directory and seed
setwd("/data/project/daylab/2019-JD-0040/")
set.seed(1234)

#Now read in the data. 
#1 = Female Saline
#2 = Female Cocaine
#3 = Male Saline
#4 = Male Cocaine
library(Seurat)
Fem_Sal_data  <- Read10X(data.dir = "1_output_Rn7_Ensembl/outs/filtered_feature_bc_matrix/")
Fem_Coc_data  <- Read10X(data.dir = "2_output_Rn7_Ensembl/outs/filtered_feature_bc_matrix/")
Male_Sal_data <- Read10X(data.dir = "3_output_Rn7_Ensembl/outs/filtered_feature_bc_matrix/")
Male_Coc_data <- Read10X(data.dir = "4_output_Rn7_Ensembl/outs/filtered_feature_bc_matrix/")

#Create the Seurat object 
#Using arbitrary cutoffs here. THis allows us to interrogate the quality of every cell, while reserving the right to remove some at 
#a later QC point. 
Fem_Sal  <- CreateSeuratObject(counts = Fem_Sal_data,min.cells = 1,min.features = 1) #6375 nuclei
Fem_Coc  <- CreateSeuratObject(counts = Fem_Coc_data,min.cells = 1,min.features = 1) #7265 nuclei
Male_Sal <- CreateSeuratObject(counts = Male_Sal_data,min.cells = 1,min.features = 1) #6213 nuclei 
Male_Coc <- CreateSeuratObject(counts = Male_Coc_data,min.cells = 1,min.features = 1) #6375 nuclei 
#26182 total nuclei

#Some of the genes in Rn7 do not contain the Mt prefix. Therefore, we are going to create a vector contain the MT_genes
#This vector will then be used to create the percent_mito feature
Rn7_gtf <- as.data.frame(rtracklayer::import("/data/project/daylab/Genome-Reference/Genomes/mRatBN7.2/Ensembl/Rattus_norvegicus.mRatBN7.2.105.gtf"))
MT_genes <- subset(Rn7_gtf,subset = (seqnames == "MT" & type == "gene"))$gene_name

#Identify the percentage of reads mapping to mitochondrial genes 
Fem_Sal <- PercentageFeatureSet(Fem_Sal, 
                                features = MT_genes[which(MT_genes %in% rownames(Fem_Sal@assays$RNA@counts))], 
                                col.name = "percent_mito",
                                assay    = "RNA")
Fem_Coc <- PercentageFeatureSet(Fem_Coc, 
                                features = MT_genes[which(MT_genes %in% rownames(Fem_Coc@assays$RNA@counts))], 
                                col.name = "percent_mito",
                                assay    = "RNA")
Male_Sal <- PercentageFeatureSet(Male_Sal, 
                                features = MT_genes[which(MT_genes %in% rownames(Male_Sal@assays$RNA@counts))], 
                                col.name = "percent_mito",
                                assay    = "RNA")
Male_Coc <- PercentageFeatureSet(Male_Coc, 
                                 features = MT_genes[which(MT_genes %in% rownames(Male_Coc@assays$RNA@counts))], 
                                 col.name = "percent_mito",
                                 assay    = "RNA")

#VlnPlots to visualize QC metrics
VlnPlot(Fem_Sal, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3,pt.size = 0)
VlnPlot(Fem_Coc, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3,pt.size = 0)
VlnPlot(Male_Sal, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3,pt.size = 0)
VlnPlot(Male_Coc, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3,pt.size = 0)

#Subset data to have greater than 200 features and less than 5% of reads mapping to mitochondrial genes 
Fem_Sal  <- subset(x = Fem_Sal, subset  =  nFeature_RNA > 200 & percent_mito < 5) #6373 nuclei 
Fem_Coc  <- subset(x = Fem_Coc, subset  =  nFeature_RNA > 200 & percent_mito < 5) #7263 nuclei 
Male_Sal <- subset(x = Male_Sal, subset =  nFeature_RNA > 200 & percent_mito < 5) #6210 nuclei 
Male_Coc <- subset(x = Male_Coc, subset =  nFeature_RNA > 200 & percent_mito < 5) #6374 nuclei
#26220 Total Nuclei

# #Replot to visualize the QC metrics following the subset 
VlnPlot(Fem_Sal, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3,pt.size = 0)
VlnPlot(Fem_Coc, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3,pt.size = 0)
VlnPlot(Male_Sal, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3,pt.size = 0)
VlnPlot(Male_Coc, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3,pt.size = 0)

#Log normalize the count data
Fem_Sal  <- NormalizeData(Fem_Sal, normalization.method = "LogNormalize", scale.factor = 10000)
Fem_Coc  <- NormalizeData(Fem_Coc, normalization.method = "LogNormalize", scale.factor = 10000)
Male_Sal <- NormalizeData(Male_Sal, normalization.method = "LogNormalize", scale.factor = 10000)
Male_Coc <- NormalizeData(Male_Coc, normalization.method = "LogNormalize", scale.factor = 10000)

#Calculat evariable features
Fem_Sal  <- FindVariableFeatures(Fem_Sal, selection.method = "vst", nfeatures = 2000)
Fem_Coc  <- FindVariableFeatures(Fem_Coc, selection.method = "vst", nfeatures = 2000)
Male_Sal <- FindVariableFeatures(Male_Sal, selection.method = "vst", nfeatures = 2000)
Male_Coc <- FindVariableFeatures(Male_Coc, selection.method = "vst", nfeatures = 2000)

#Add metadata characteristics
Fem_Sal$Stim  <- "Saline"
Fem_Coc$Stim  <- "Cocaine"
Male_Sal$Stim <- "Saline"
Male_Coc$Stim <- "Cocaine"

Fem_Sal$Sex   <- "Female"
Fem_Coc$Sex   <- "Female"
Male_Sal$Sex  <- "Male"
Male_Coc$Sex  <- "Male"

Fem_Sal$Sex_Stim  <- "Fem_Sal"
Fem_Coc$Sex_Stim  <- "Fem_Coc"
Male_Sal$Sex_Stim <- "Male_Sal"
Male_Coc$Sex_Stim <- "Male_Coc"

Fem_Sal$GEM  <- 1
Fem_Coc$GEM  <- 2
Male_Sal$GEM <- 3
Male_Coc$GEM <- 4

Fem_Sal$Dataset  <- 1
Fem_Coc$Dataset  <- 2
Male_Sal$Dataset <- 3
Male_Coc$Dataset <- 4

#Integrate all datasets 
#Will use 17 dimensions and 0.2 resolution. 
All_Groups_log <- FindIntegrationAnchors(object.list = list(Fem_Sal,Fem_Coc,Male_Sal,Male_Coc), dims = 1:17)
All_Groups_log  <- IntegrateData(anchorset = All_Groups_log,dims = 1:17)

DefaultAssay(All_Groups_log) <- "integrated"

# Run the standard workflow for visualization and clustering
All_Groups_log <- ScaleData(All_Groups_log,verbose = FALSE)
All_Groups_log <- RunPCA(All_Groups_log,npcs = 17,verbose = FALSE) #Compute 50 npcs by default
# Dimensionality reduction and Clustering
All_Groups_log <- RunUMAP(All_Groups_log, reduction = "pca", dims = 1:17)
All_Groups_log <- FindNeighbors(All_Groups_log, reduction = "pca", dims = 1:17)

All_Groups_log <- FindClusters(All_Groups_log, resolution = 0.1)

#Make the umap
DimPlot(All_Groups_log,label = TRUE) + NoLegend()

#Save the umap as v1
library(ggplot2)
ggsave(plot = DimPlot(All_Groups_log,label = TRUE) + NoLegend(),
       filename = "MCN_Code/Plots/Repeated_Cocaine_umap_v1.pdf")


#Now identify the cell types
Dp1 <- DotPlot(object = All_Groups_log,
               assay  = "RNA",
               features =  c("Drd1","Pdyn","Ebf1", #D1-MSN
                             "Drd2","Penk", #D2-MSN
                             "Drd3","Grm8", #D3/Grm8 MSN
                             "Elavl2","Kit","Sst", #GABAergic markers
                             "Slc17a7", #Glut
                             "Mbp","Opalin", #Oligs
                             "Aqp4","Gja1", #Astrocytes
                             "Pdgfra", #Polydendrocytes
                             "Arhgap15", #Microglia
                             "Rgs5", #Mural
                             "Ppp1r1b","Foxp2","Bcl11b","Gad1","Syt1"), #Neuronal markers
               cols = c("lightgrey","red")) +
  theme(axis.text.x = element_text(angle = 45, hjust =1))

ggsave(plot     = Dp1,
       filename = "MCN_Code/Plots/DotPlot_DoubletsMaybe_v1.pdf",
       height   = 12,
       width    = 12)


#Feature plot
DefaultAssay(All_Groups_log) <- "RNA"
Fp1 <- FeaturePlot(object = All_Groups_log,
            features =  c("Drd1","Pdyn","Ebf1", #D1-MSN
                          "Drd2","Penk", #D2-MSN
                          "Drd3","Grm8", #D3/Grm8 MSN
                          "Elavl2","Kit","Sst", #GABAergic markers
                          "Slc17a7", #Glut
                          "Mbp","Opalin", #Oligs
                          "Aqp4","Gja1", #Astrocytes
                          "Pdgfra", #Polydendrocytes
                          "Arhgap15", #Microglia
                          "Rgs5", #Mural
                          "Ppp1r1b","Foxp2","Bcl11b","Gad1","Syt1"), #Neuronal markers
            cols = c("lightgrey","red"))

ggsave(plot     = Fp1,
       filename = "MCN_Code/Plots/FeaturePlots_v1.pdf",
       height   = 20,
       width    = 20)

# 0	Olig-1
# 1	Drd1-MSN
# 2	Drd2-MSN
# 3	GABAergic-1
# 4	Grm8-MSN
# 5	Astrocytes
# 6	Polydendrocyte
# 7	Microglia
# 8	Glutamatergic
# 9	Drd2-MSN-2
# 10	Drd3-MSN
# 11	Sst-Interneuron
# 12	Doublets
# 13	Mural
# 14	Ebf1-pos
# 15	GABAergic-2

#Rename the identities
All_Groups_log <-  RenameIdents(object = All_Groups_log,
                                "0" = "Olig-1",#
                                "1" = "Drd1-MSN",#
                                "2" = "Drd2-MSN",#
                                "3" = "GABAergic-1",#
                                "4" = "Grm8-MSN",#
                                "5" = "Astrocyte",#
                                "6" = "Polydendrocyte",#
                                "7" = "Microglia",
                                "8" = "Glutamatergic",#
                                "9" = "Drd2-MSN-2",#
                                "10" = "Drd3-MSN",#
                                "11" = "Sst-Interneuron",#
                                "12" = "Doublets",#
                                "13" = "Mural",#
                                "14" = "Pos-Ebf1",#
                                "15" = "GABAergic-2")

#Create a cell type column
All_Groups_log$CellType <- Idents(All_Groups_log)


DimPlot(object = All_Groups_log,reduction = "umap",label = TRUE) + 
  NoLegend()

ggsave(plot = DimPlot(All_Groups_log,label = TRUE) + NoLegend(),
       filename = "MCN_Code/Plots/Repeated_Cocaine_umap_v1_withCellTypeLabels.pdf")


#save the object
saveRDS(object = All_Groups_log,file = "MCN_Code/Objects/Repeated_Cocaine_v1.RDS")

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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_3.3.2      SeuratObject_4.0.2 Seurat_4.0.4      
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15                  colorspace_1.4-1            deldir_1.0-6                ellipsis_0.3.2             
# [5] ggridges_0.5.2              XVector_0.28.0              GenomicRanges_1.40.0        rstudioapi_0.11            
# [9] spatstat.data_2.1-0         farver_2.0.3                leiden_0.3.3                listenv_0.8.0              
# [13] ggrepel_0.8.2               RSpectra_0.16-0             codetools_0.2-16            splines_4.0.2              
# [17] polyclip_1.10-0             jsonlite_1.7.1              Rsamtools_2.4.0             ica_1.0-2                  
# [21] cluster_2.1.0               png_0.1-7                   uwot_0.1.10                 shiny_1.5.0                
# [25] sctransform_0.3.2           spatstat.sparse_2.0-0       compiler_4.0.2              httr_1.4.2                 
# [29] Matrix_1.3-4                fastmap_1.1.0               lazyeval_0.2.2              cli_3.3.0                  
# [33] later_1.1.0.1               htmltools_0.5.2             tools_4.0.2                 igraph_1.2.6               
# [37] gtable_0.3.0                glue_1.6.2                  GenomeInfoDbData_1.2.3      RANN_2.6.1                 
# [41] reshape2_1.4.4              dplyr_1.0.2                 rappdirs_0.3.1              Rcpp_1.0.7                 
# [45] Biobase_2.48.0              scattermore_0.7             vctrs_0.4.1                 Biostrings_2.56.0          
# [49] nlme_3.1-149                rtracklayer_1.48.0          lmtest_0.9-38               stringr_1.4.0              
# [53] globals_0.13.1              mime_0.9                    miniUI_0.1.1.1              lifecycle_0.2.0            
# [57] irlba_2.3.3                 XML_3.99-0.5                goftest_1.2-2               future_1.19.1              
# [61] MASS_7.3-53                 zlibbioc_1.34.0             zoo_1.8-8                   scales_1.1.1               
# [65] spatstat.core_2.3-0         promises_1.1.1              spatstat.utils_2.2-0        SummarizedExperiment_1.18.2
# [69] parallel_4.0.2              RColorBrewer_1.1-2          yaml_2.2.1                  reticulate_1.16            
# [73] pbapply_1.4-3               gridExtra_2.3               rpart_4.1-15                stringi_1.7.6              
# [77] S4Vectors_0.26.1            BiocGenerics_0.34.0         BiocParallel_1.22.0         GenomeInfoDb_1.24.2        
# [81] rlang_1.0.2                 pkgconfig_2.0.3             matrixStats_0.57.0          bitops_1.0-6               
# [85] lattice_0.20-41             ROCR_1.0-11                 purrr_0.3.4                 tensor_1.5                 
# [89] labeling_0.3                GenomicAlignments_1.24.0    patchwork_1.0.1             htmlwidgets_1.5.2          
# [93] cowplot_1.1.0               tidyselect_1.1.0            RcppAnnoy_0.0.19            plyr_1.8.6                 
# [97] magrittr_1.5                R6_2.4.1                    IRanges_2.22.2              generics_0.0.2             
# [101] DelayedArray_0.14.1         withr_2.3.0                 pillar_1.4.6                mgcv_1.8-33                
# [105] fitdistrplus_1.1-1          survival_3.2-7              abind_1.4-5                 RCurl_1.98-1.2             
# [109] tibble_3.0.4                future.apply_1.6.0          crayon_1.3.4                KernSmooth_2.23-17         
# [113] spatstat.geom_2.4-0         plotly_4.9.2.1              grid_4.0.2                  data.table_1.13.0          
# [117] digest_0.6.26               xtable_1.8-4                tidyr_1.1.2                 httpuv_1.5.4               
# [121] stats4_4.0.2                munsell_0.5.0               viridisLite_0.3.0   
