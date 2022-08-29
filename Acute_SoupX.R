#Now that I have built a general acute cocaine object. I am going to run SoupX. 
#If Soupx doesn't fix any problems, I will then run doublet finder. 

#First set working directory and seed. 
setwd("/data/project/daylab/2019-JD-0037/")
set.seed(1234)

#load seurat and soupx. 
library(Seurat)
library(SoupX)
library(ggplot2)

#Load the object
All_Groups_log <- readRDS(file = "/data/project/daylab/2019-JD-0040/MCN_Code/Objects/Acute_Cocaine_v1.RDS")

#Check the map
DimPlot(All_Groups_log,reduction = ,label = TRUE) + NoLegend()

#Create a cell column
All_Groups_log$Cell <- row.names(All_Groups_log@meta.data)

#Create a non-specific cell id
All_Groups_log$Non_Specific_CellID <- as.character(lapply(strsplit(x = All_Groups_log$Cell,split = "_"),"[",1))

#Create an empty vector that will hold the matrices and single GEM objects
All_Groups_adjusted_matrices <- vector(mode = "list",length = 4)
All_Groups_single_objects    <- vector(mode = "list",length = 4)

#Get all of the directories holding cellranger output. 
Data_Dirs <- factor(x = list.files(path = ".")[grep("output_Rn7_Ensembl",list.files(path = "."))],
                    levels = c("FemSal_output_Rn7_Ensembl",
                               "FemCoc_output_Rn7_Ensembl",
                               "MaleSal_output_Rn7_Ensembl",
                               "MaleCoc_output_Rn7_Ensembl"))



#Loop through to run the standard SoupX workflow
for(i in 1:4){
  #Read counts and droplet matrices
  Counts <- Seurat::Read10X(file.path(paste0(levels(Data_Dirs)[i],"/outs/filtered_feature_bc_matrix/")))
  Droplets <- Seurat::Read10X(file.path(paste0(levels(Data_Dirs)[i],"/outs/raw_feature_bc_matrix/")))
  #Create soup channel
  sc  <- SoupChannel(tod = Droplets, toc = Counts)
  #Add Metadata info that includes cell and celltype
  sc$metaData$Cell <- row.names(sc$metaData)
  print("Starting to add metadata")
  start <- Sys.time()
  print(start)
  sc$metaData$CellType <- apply(X = sc$metaData,MARGIN = 1,function(x){
    as.character(subset(All_Groups_log@meta.data,
                        subset=(Non_Specific_CellID == x["Cell"] & GEM == i))$CellType)
  })
  end <- Sys.time()
  print(end)
  print(end-start)
  print("Done adding metadata")
  print("Starting to set clusters")
  #Finalize the SoupX workflow
  sc <- setClusters(sc,clusters = as.character(sc$metaData$CellType))
  sc  <- autoEstCont(sc)
  out <- adjustCounts(sc,roundToInt = TRUE)
  All_Groups_adjusted_matrices[[i]] <- out
  meta <-  All_Groups_log@meta.data[which(All_Groups_log$Non_Specific_CellID %in% colnames(All_Groups_adjusted_matrices[[i]]) & All_Groups_log$GEM == i),]
  row.names(meta) <- meta$Non_Specific_CellID
  All_Groups_single_objects[[i]] <- CreateSeuratObject(All_Groups_adjusted_matrices[[i]],
                                                       meta.data = meta)
  print(paste(i,"is finished"))
  rm(sc,out,Counts,Droplets)
}


#Some of the genes in Rn7 do not contain the Mt prefix. Therefore, we are going to create a vector contain the MT_genes
#This vector will then be used to create the percent_mito feature
Rn7_gtf <- as.data.frame(rtracklayer::import("/data/project/daylab/Genome-Reference/Genomes/mRatBN7.2/Ensembl/Rattus_norvegicus.mRatBN7.2.105.gtf"))
MT_genes <- subset(Rn7_gtf,subset = (seqnames == "MT" & type == "gene"))$gene_name


All_Groups_log <- All_Groups_single_objects

#Now calculate the mitochondrial percentage for every cell
All_Groups_log <- lapply(X   = All_Groups_log,
                         FUN = function(x){
                           PercentageFeatureSet(x, 
                                                features = MT_genes[which(MT_genes %in% rownames(x@assays$RNA@counts))], 
                                                col.name = "percent_mito",
                                                assay    = "RNA")
                         })


#Plot the number of genes, UMIs, and percent_mito for every cell
lapply(X   = All_Groups_log,
       FUN = function(x){
         VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3,pt.size = 0)
       })


#Now subset for cells with > 200 genes and less than 5% of reads mapping to the mitochondrial genome
All_Groups_log <- lapply(X   = All_Groups_log,
                         FUN = function(x){
                           subset(x = x, subset  =  nFeature_RNA > 200 & percent_mito < 5)
                         })

#Perform log normalization
All_Groups_log <- lapply(X   = All_Groups_log,
                         FUN = function(x){
                           NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
                         })

#Find variable features
All_Groups_log <- lapply(X   = All_Groups_log,
                         FUN = function(x){
                           FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
                         })


#Integrate the data
All_Groups_log_Decontaminated  <- FindIntegrationAnchors(object.list = All_Groups_log, dims = 1:17)
All_Groups_log_Decontaminated  <- IntegrateData(anchorset = All_Groups_log_Decontaminated,dims = 1:17)

# Run the standard workflow for visualization and clustering
All_Groups_log_Decontaminated <- ScaleData(All_Groups_log_Decontaminated,verbose = FALSE)
# Dimensionality reduction and Clustering
All_Groups_log_Decontaminated <- RunPCA(All_Groups_log_Decontaminated,npcs = 17 ,verbose = FALSE) #Compute 50 npcs by default
All_Groups_log_Decontaminated <- RunUMAP(All_Groups_log_Decontaminated, reduction = "pca", dims = 1:17)
All_Groups_log_Decontaminated <- FindNeighbors(All_Groups_log_Decontaminated, reduction = "pca", dims = 1:17)

#Find clusters
#DefaultAssay(VTA_Decontaminated) <- "integrated"
All_Groups_log_Decontaminated <- FindClusters(All_Groups_log_Decontaminated, resolution = 0.1)


DimPlot(All_Groups_log_Decontaminated,label = TRUE,reduction = "umap") + NoLegend()

library(ggplot2)
Dp1 <- DotPlot(object = All_Groups_log_Decontaminated,
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
       filename = "/data/project/daylab/2019-JD-0040/MCN_Code/Plots/Acute_Doplot_v2_postSoupX.pdf",
       height   = 12,
       width    = 12)


#Feature plot
DefaultAssay(All_Groups_log_Decontaminated) <- "RNA"
Fp1 <- FeaturePlot(object = All_Groups_log_Decontaminated,
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
       filename = "/data/project/daylab/2019-JD-0040/MCN_Code/Plots/FeaturePlots_Acute_v2_postSoupX.pdf",
       height   = 20,
       width    = 20)


#0 Olig-1
#1 Drd1-MSN
#2 Drd2-MSN
#3 GABAergic-1
#4 Astrocyte
#5 Grm8-MSN
#6 Polydendrocyte
#7 Microglia
#8 Pvalb-Interneuron
#9 Drd3-MSN
#10 Sst-Interneuron
#11 Glutamatergic
#12 Mural

#Rename the identities
All_Groups_log_Decontaminated <-  RenameIdents(object = All_Groups_log_Decontaminated,
                                               "0" = "Olig-1",#
                                               "1" = "Drd1-MSN",#
                                               "2" = "Drd2-MSN",#
                                               "3" = "GABAergic-1",#
                                               "4" = "Astrocyte",#
                                               "5" = "Grm8-MSN",#
                                               "6" = "Polydendrocyte",#
                                               "7" = "Microglia",
                                               "8" = "Pvalb-Interneuron",#
                                               "9" = "Drd3-MSN",#
                                               "10" = "Sst-Interneuron",#
                                               "11" = "Glutamatergic",#
                                               "12" = "Mural")

D_Umap <- DimPlot(All_Groups_log_Decontaminated,label = TRUE,reduction = "umap") + 
  NoLegend() +
  ggtitle("Decontaminated") +
  theme(plot.title = element_text(hjust = 0.5))


All_Groups_log <- readRDS(file = "/data/project/daylab/2019-JD-0040/MCN_Code/Objects/Acute_Cocaine_v1.RDS")
O_Umap <- DimPlot(All_Groups_log,label = TRUE,reduction = "umap") + 
  NoLegend() +
  ggtitle("Original") +
  theme(plot.title = element_text(hjust = 0.5))
  
ggsave(plot = cowplot::plot_grid(plotlist = list(O_Umap,D_Umap),nrow = 1),
       filename = "/data/project/daylab/2019-JD-0040/MCN_Code/Plots/Orig_Decont_Acute_Umaps.pdf",
       height = 12,
       width = 12)

saveRDS(object = All_Groups_log_Decontaminated,
        file = "/data/project/daylab/2019-JD-0040/MCN_Code/Objects/Acute_Decontaminated.RDS")


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
#   [1] ggplot2_3.3.2      SoupX_1.5.2        SeuratObject_4.0.2 Seurat_4.0.4      
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15                  colorspace_1.4-1            deldir_1.0-6                ellipsis_0.3.2             
# [5] ggridges_0.5.2              XVector_0.28.0              GenomicRanges_1.40.0        rstudioapi_0.11            
# [9] spatstat.data_2.1-0         leiden_0.3.3                listenv_0.8.0               farver_2.0.3               
# [13] ggrepel_0.8.2               RSpectra_0.16-0             codetools_0.2-16            splines_4.0.2              
# [17] polyclip_1.10-0             jsonlite_1.7.1              Rsamtools_2.4.0             ica_1.0-2                  
# [21] cluster_2.1.0               png_0.1-7                   uwot_0.1.10                 shiny_1.5.0                
# [25] sctransform_0.3.2           spatstat.sparse_2.0-0       compiler_4.0.2              httr_1.4.2                 
# [29] Matrix_1.3-4                fastmap_1.1.0               lazyeval_0.2.2              cli_3.3.0                  
# [33] later_1.1.0.1               htmltools_0.5.2             tools_4.0.2                 igraph_1.2.6               
# [37] GenomeInfoDbData_1.2.3      gtable_0.3.0                glue_1.6.2                  RANN_2.6.1                 
# [41] reshape2_1.4.4              dplyr_1.0.2                 rappdirs_0.3.1              Rcpp_1.0.7                 
# [45] Biobase_2.48.0              scattermore_0.7             Biostrings_2.56.0           vctrs_0.4.1                
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
# [89] GenomicAlignments_1.24.0    patchwork_1.0.1             htmlwidgets_1.5.2           labeling_0.3               
# [93] cowplot_1.1.0               tidyselect_1.1.0            RcppAnnoy_0.0.19            plyr_1.8.6                 
# [97] magrittr_1.5                R6_2.4.1                    IRanges_2.22.2              generics_0.0.2             
# [101] DelayedArray_0.14.1         pillar_1.4.6                withr_2.3.0                 mgcv_1.8-33                
# [105] fitdistrplus_1.1-1          survival_3.2-7              abind_1.4-5                 RCurl_1.98-1.2             
# [109] tibble_3.0.4                future.apply_1.6.0          crayon_1.3.4                KernSmooth_2.23-17         
# [113] spatstat.geom_2.4-0         plotly_4.9.2.1              grid_4.0.2                  data.table_1.13.0          
# [117] digest_0.6.26               xtable_1.8-4                tidyr_1.1.2                 httpuv_1.5.4               
# [121] stats4_4.0.2                munsell_0.5.0               viridisLite_0.3.0