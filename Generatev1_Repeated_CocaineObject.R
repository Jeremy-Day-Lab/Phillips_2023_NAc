#Data from this analysis comes from male and female adult Sprague-Dawley rats that were treated for 7 consecutive days 
#with 20mg/kg of cocaine. 1 hour forlloiwng final exposure, animals were sacrificed, and NAc was taken. 
#Goal of this analysis: Generate a seurat object for this dataset. We will later use SoupX to remove transcripts found overrepresented
#in empty droplets, and doubletfinder to remove doublets. 

#Set working directory and seed
setwd("/data/project/daylab/2019-JD-0040")
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

#Identify the percentage of reads mapping to mitochondrial genes 
Fem_Sal  <- PercentageFeatureSet(Fem_Sal, pattern = "^Mt-", col.name = "percent_mito")
Fem_Coc  <- PercentageFeatureSet(Fem_Coc, pattern = "^Mt-", col.name = "percent_mito")
Male_Sal <- PercentageFeatureSet(Male_Sal, pattern = "^Mt-", col.name = "percent_mito")
Male_Coc <- PercentageFeatureSet(Male_Coc, pattern = "^Mt-", col.name = "percent_mito")


#VlnPlots to visualize QC metrics
VlnPlot(Fem_Sal, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3,pt.size = 0)
VlnPlot(Fem_Coc, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3,pt.size = 0)
VlnPlot(Male_Sal, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3,pt.size = 0)
VlnPlot(Male_Coc, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3,pt.size = 0)

#Subset data to have greater than 200 features and less than 5% of reads mapping to mitochondrial genes 
Fem_Sal  <- subset(x = Fem_Sal, subset  =  nFeature_RNA > 200 & percent_mito < 5) #6373 nuclei 
Fem_Coc  <- subset(x = Fem_Coc, subset  =  nFeature_RNA > 200 & percent_mito < 5) #7263 nuclei 
Male_Sal <- subset(x = Male_Sal, subset =  nFeature_RNA > 200 & percent_mito < 5) #6213 nuclei 
Male_Coc <- subset(x = Male_Coc, subset =  nFeature_RNA > 200 & percent_mito < 5) #6375 nuclei
#26175 Total Nuclei

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

All_Groups_log <- FindClusters(All_Groups_log, resolution = 0.2)

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
# 3	Grm8-MSN
# 4	GABAergic-1
# 5	Astrocye
# 6	Polydendrocyte
# 7	Microglia
# 8	Pvalb-Interneuron
# 9	Glutamatergic
# 10	Drd2-MSN-2
# 11	Drd3-MSN
# 12	Doublets
# 13	Sst-Interneuron
# 14	Doublets-2
# 15	Olig-2
# 16	Mural
# 17	Ebf1-pos
# 18	GABAergic-2

#Rename the identities
All_Groups_log <-  RenameIdents(object = All_Groups_log,
                                "0" = "Olig-1",#
                                "1" = "Drd1-MSN",#
                                "2" = "Drd2-MSN",#
                                "3" = "Grm8-MSN",#
                                "4" = "GABAergic-1",#
                                "5" = "Astrocyte",#
                                "6" = "Polydendrocyte",#
                                "7" = "Microglia",
                                "8" = "Pvalb-Interneuron",#
                                "9" = "Glutamatergic",#
                                "10" = "Drd2-MSN-2",#
                                "11" = "Drd3-MSN",#
                                "12" = "Doublets",#
                                "13" = "Sst-Interneuron",#
                                "14" = "Doublets-2",#
                                "15" = "Olig-2",
                                "16" = "Mural",
                                "17" = "Ebf1-Pos",
                                "18" = "GABAergic-2")

#Create a cell type column
All_Groups_log$CellType <- Idents(All_Groups_log)


DimPlot(object = All_Groups_log,reduction = "umap",label = TRUE) + 
  NoLegend()

#
