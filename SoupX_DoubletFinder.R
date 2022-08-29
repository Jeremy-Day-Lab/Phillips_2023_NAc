#Now that I have built a general repeated cocaine object. I am going to run SoupX. 
#If Soupx doesn't fix any problems, I will then run doublet finder. 

#First set working directory and seed. 
setwd("/data/project/daylab/2019-JD-0040")
set.seed(1234)

#load seurat and soupx. 
library(Seurat)
library(SoupX)
library(ggplot2)

#Load the object
All_Groups_log <- readRDS(file = "MCN_Code/Objects/Repeated_Cocaine_v1.RDS")

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
Data_Dirs <- list.files(path = ".")[grep("output_Rn7_Ensembl",list.files(path = "."))]

#Loop through to run the standard SoupX workflow
for(i in 1:4){
  print(Data_Dirs[i])
  #Read counts and droplet matrices
  Counts <- Seurat::Read10X(file.path(paste0(Data_Dirs[i],"/outs/filtered_feature_bc_matrix/")))
  Droplets <- Seurat::Read10X(file.path(paste0(Data_Dirs[i],"/outs/raw_feature_bc_matrix/")))
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
  out <- adjustCounts(sc)
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

#Build a dotplot for cluster characterization
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
       filename = "MCN_Code/Plots/DotPlot_PostDecontamination_v1.pdf",
       height   = 12,
       width    = 12)

#Rename the identities
All_Groups_log_Decontaminated <-  RenameIdents(object = All_Groups_log_Decontaminated,
                                               "0" = "Olig-1",#
                                               "1" = "Drd1-MSN",#
                                               "2" = "Drd2-MSN",#
                                               "3" = "GABAergic-1",#
                                               "4" = "Grm8-MSN",#
                                               "5" = "Astrocyte",#
                                               "6" = "Polydendrocyte",#
                                               "7" = "Microglia",#
                                               "8" = "Glutamatergic",#
                                               "9" = "Drd2-MSN-2",#
                                               "10" = "Drd3-MSN",#
                                               "11" = "Sst-Interneuron",#
                                               "12" = "Doublets",#
                                               "13" = "Mural",#
                                               "14" = "Pos-Ebf1")#

#Relevel the identities
All_Groups_log_Decontaminated$CellType <- factor(x = All_Groups_log_Decontaminated$CellType,
                                                levels = c("Drd1-MSN",#
                                                           "Drd2-MSN",#
                                                           "Drd2-MSN-2",#
                                                           "Drd3-MSN",#
                                                           "GABAergic-1",#
                                                           "Grm8-MSN",#
                                                           "Glutamatergic",#
                                                           "Sst-Interneuron",#
                                                           "Astrocyte",#
                                                           "Microglia",#
                                                           "Mural",#
                                                           "Olig-1",#
                                                           "Polydendrocyte",#
                                                           "Pos-Ebf1",#
                                                           "Doublets"))

Idents(All_Groups_log_Decontaminated) <- All_Groups_log_Decontaminated$CellType

#Read in the original object and relevel identities
All_Groups_log_original <- readRDS(file = "MCN_Code/Objects/Repeated_Cocaine_v1.RDS")

All_Groups_log_original$CellType <- factor(x = All_Groups_log_original$CellType,
                                           levels = c("Drd1-MSN",#
                                                      "Drd2-MSN",#
                                                      "Drd2-MSN-2",#
                                                      "Drd3-MSN",#
                                                      "GABAergic-1",#
                                                      "GABAergic-2",
                                                      "Grm8-MSN",#
                                                      "Glutamatergic",#
                                                      "Sst-Interneuron",#
                                                      "Astrocyte",#
                                                      "Microglia",#
                                                      "Mural",#
                                                      "Olig-1",
                                                      "Polydendrocyte",#
                                                      "Pos-Ebf1",#
                                                      "Doublets"))

Idents(All_Groups_log_original) <- All_Groups_log_original$CellType


#Compare Oligodendrocyte genes between the two objects. These tend to be problematic
VP_O <- VlnPlot(object = All_Groups_log_original,features = c("Mbp","Mobp","Opalin"),assay = "RNA",stack = TRUE,flip = TRUE) + 
  NoLegend() + 
  ggtitle("Original") +
  theme(plot.title = element_text(hjust = 0.5))
VP_D <- VlnPlot(object = All_Groups_log_Decontaminated,features = c("Mbp","Mobp","Opalin"),assay = "RNA",stack = TRUE,flip = TRUE) + 
  NoLegend() + 
  ggtitle("Decontaminated") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = cowplot::plot_grid(plotlist = list(VP_O,VP_D),nrow = 1),
       filename = "./MCN_Code/Plots/VP_Comp_Orig_vs_Decontaminated_OligGenes.pdf")

library(ggplot2)
O_Umap <- DimPlot(All_Groups_log_original,label = TRUE,reduction = "umap") + 
  NoLegend() + 
  ggtitle("Original Object") +
  theme(plot.title = element_text(hjust = 0.5))

D_Umap <- DimPlot(All_Groups_log_Decontaminated,label = TRUE,reduction = "umap") + 
  NoLegend() + 
  ggtitle("Decontaminated Object") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = cowplot::plot_grid(plotlist = list(O_Umap,D_Umap),nrow = 1),
       filename = "./MCN_Code/Plots/Umap_Comp_Orig_vs_Decontaminated.pdf")


saveRDS(object = All_Groups_log_Decontaminated,file = "MCN_Code/Objects/Repeated_PostSoupX.RDS")


######DOUBLET FINDER#######
#Still seem to be doublets identified in these objects. Going to run doublet finder
library(DoubletFinder)

#Change the name of the object to All_Groups_Log_Repeated 
All_Groups_Log_Repeated <- All_Groups_log_Decontaminated

#Now that every cell has a celltype, I need to go back into each individual object to identify doublets. 
#Doublet finder needs the cells to have a cell type because the final number of detectable doublets by multiplying the proportion of homotypic doublets 
#(which is defined as the sum of squared cell state frequencies) by the total doublet rate, which I found in the chromium chemistry v3 handbook. 
#pN = Number of artificial doublets used during dimensionality reduciton and identification of doublets.
#pK = Neighborhood size 

#Make a column in the All_Groups_metadata that does contains non-unique cell names. Seruat forces unique cell names by adding an underscore and then an integer. 
#Thus, splitting by the underscore gives us the original cell identity that can be matched back to the original Seurat object. 
All_Groups_Log_Repeated$Non_unique_Cell_Name <- as.character(lapply(strsplit(rownames(All_Groups_Log_Repeated@meta.data),split = "_"),"[",1))

#Make a column that identifies the forced unique cell names by writing out the rownames
All_Groups_Log_Repeated$Unique_Cell_Name <- row.names(All_Groups_Log_Repeated@meta.data)

##########FEMALE SALINE##########
#Pull the FeMale_Sal_Repeatedine metadata from the main object. Here I only take the Non unique cell name created in the first command, the cell type, and the unique
#cell name forced by Seurat
FS_metadata <- subset(All_Groups_Log_Repeated,subset = (Sex_Stim == "Fem_Sal"))@meta.data[,c("Non_unique_Cell_Name","CellType","Unique_Cell_Name")]


#Create Fem_Sal_Repeated from the list of single objects generated above. 
Fem_Sal_Repeated <- All_Groups_log[[1]]


#Fem_Sal_Repeated should have rownames that correspond to the non-unique_cell names. Now, I make that a column which will be used to merge in the following command 
Fem_Sal_Repeated$Non_unique_Cell_Name <- rownames(Fem_Sal_Repeated@meta.data)

#Merge the two dataframes by non-unique-cellname
FS <- merge(x  = Fem_Sal_Repeated@meta.data,
            y  = FS_metadata,
            by = "Non_unique_Cell_Name")

#Make the rownames the non_unique cellname
rownames(FS) <- FS$Non_unique_Cell_Name

#Make sure nothing is duplicated when merging. If nothing is duplicated then the number of unique rows should be the same as the number of rows when pulling Female Saline 
#cells from the All_Groups_Log_Repeated metadata
library(dplyr)
nrow(FS %>% distinct()) == nrow(FS_metadata) #TRUE

#Now I want make sure FS to be in the same order as Fem_Sal_Repeated. I will use the match function as demonstrated above 
FS <- FS[match(row.names(Fem_Sal_Repeated@meta.data),row.names(FS)),]
#Double check that they are in the same order, and that they are the same length
all(row.names(FS) == row.names(Fem_Sal_Repeated@meta.data)) #TRUE
all(length(row.names(FS)) == length(row.names(Fem_Sal_Repeated@meta.data))) #TRUE

# #Add the cell type to the metadata
# Fem_Sal_Repeated <- AddMetaData(object = Fem_Sal_Repeated,metadata = FS$CellType,col.name = "CellType")

#DoubletFinder requires that each individual object is scaled, and undergoes standard dimensionality reduction. 
#For this, I use the same number of PCs used in the processing of the integrated object.Fem_Sal_Repeated <- ScaleData(Fem_Sal_Repeated,verbose = FALSE)
Fem_Sal_Repeated <- ScaleData(Fem_Sal_Repeated,verbose = FALSE)
Fem_Sal_Repeated <- RunPCA(Fem_Sal_Repeated,npcs = 17,verbose = FALSE)
Fem_Sal_Repeated <- RunUMAP(Fem_Sal_Repeated,reduction = "pca",dims = 1:17)

#paramSweep_v3 tests multiple pN and pK parameters to identify the optimal pK value. The Cell Systems paper the authors of Doublet finder found that 
#pK is "the main parameter that needs to be tuned", and that results are largely invariate of the pN value. 
Fem_Sal_Repeated_sweep <- paramSweep_v3(Fem_Sal_Repeated,PCs = 1:17,sct = FALSE)

#Summarizes the sweep. GT stands for ground truth which we do not have here. 
Fem_Sal_Repeated_sweep_stats <- summarizeSweep(sweep.list = Fem_Sal_Repeated_sweep,GT = FALSE)

#Identify and plot the pK 
FS_pk <- find.pK(sweep.stats = Fem_Sal_Repeated_sweep_stats)

#Identify the pK value corresponding to the maximum mean-variance-normalized-bimodality coefficient
FS_pk <- as.numeric(as.character(FS_pk[which(FS_pk$BCmetric == max(FS_pk$BCmetric)),"pK"]))

#Now model the number of homotypic doublets. This calculates the sum of the squares of the proportions of the annotations. 
FS_homotypic_props <- modelHomotypic(Fem_Sal_Repeated$CellType)

#The first calculation multiplies the number of cells by the doublet rate
nExp_poi <- round(.10*nrow(Fem_Sal_Repeated@meta.data))
#Now the number of expected homotypic doublets is multiplied by the proportion of Heterotypic doublets 
nExp_poi_adj <- round(nExp_poi*(1-FS_homotypic_props))

#Find Doublets using the less strict nExp_poi_adj. This is because DoubletFinder cannot confidently identify homotypic doublets, and using nExp_poi may 
#remove singlets that the program identifies as doublets
#Fem_Sal_Repeated <- doubletFinder_v3(seu = Fem_Sal_Repeated,PCs = 1:17,pN = 0.25,pK = FS_pk,nExp = nExp_poi_adj,reuse.pANN = FALSE,sct = FALSE)
Fem_Sal_Repeated <- doubletFinder_v3(seu = Fem_Sal_Repeated,PCs = 1:17,pN = 0.25,pK = FS_pk,nExp = nExp_poi,reuse.pANN = FALSE,sct = FALSE)
#change names of the dataframe so everything can be merged later
names(Fem_Sal_Repeated@meta.data)[c(16,17)] <- c("pANN","DF.Classification")

#remove all objects except for Fem_Sal_Repeated so that environment does not become overloaded 
rm(FS,FS_metadata,FS_homotypic_props,FS_pk,nExp_poi,nExp_poi_adj,Fem_Sal_Repeated_sweep,Fem_Sal_Repeated_sweep_stats)

##########FEMALE Cocaine##########
#Moving forward the same workflow will be used but extensive explanation will not be provided. 
#Pull the FeMale_Coc_Repeatedine metadata from the main object that contains the cellnames 
FC_metadata <- subset(All_Groups_Log_Repeated,subset = (Sex_Stim == "Fem_Coc"))@meta.data[,c("Non_unique_Cell_Name","CellType","Unique_Cell_Name")]

#Create Fem_Sal_Repeated from the list of single objects generated above. 
Fem_Coc_Repeated <- All_Groups_log[[2]]

#Make the rownames a column 
Fem_Coc_Repeated$Non_unique_Cell_Name <- rownames(Fem_Coc_Repeated@meta.data)

#Merge the two dataframes by non-unique-cellname
FC <- merge(x  = Fem_Coc_Repeated@meta.data,
            y  = FC_metadata,
            by = "Non_unique_Cell_Name")

#Make non unique names the rownames
rownames(FC) <- FC$Non_unique_Cell_Name

#Make sure nothing is duplicated 
nrow(FC %>% distinct()) == nrow(FC_metadata) #TRUE 

#Make sure these are in the same order
#I want FC to be in the same order as Fem_Coc_Repeated, therefore FC comes second and Fem_Coc_Repeated comes first 
FC<- FC[match(row.names(Fem_Coc_Repeated@meta.data),row.names(FC)),]

#Check everything matches
all(row.names(FC) == row.names(Fem_Coc_Repeated@meta.data))
all(length(row.names(FC)) == length(row.names(Fem_Coc_Repeated@meta.data)))

# #Add the metadata 
# Fem_Coc_Repeated <- AddMetaData(object = Fem_Coc_Repeated,metadata = FC$CellType,col.name = "CellType")

#Scale and dimensionality reduction
Fem_Coc_Repeated <- ScaleData(Fem_Coc_Repeated,verbose = FALSE)
Fem_Coc_Repeated <- RunPCA(Fem_Coc_Repeated,npcs = 17,verbose = FALSE)
Fem_Coc_Repeated <- RunUMAP(Fem_Coc_Repeated,reduction = "pca",dims = 1:17)

#Sweep PNs and PKs
Fem_Coc_Repeated_sweep <- paramSweep_v3(Fem_Coc_Repeated,PCs = 1:17,sct = FALSE)

#Summarize the sweekp
Fem_Coc_Repeated_sweep_stats <- summarizeSweep(sweep.list = Fem_Coc_Repeated_sweep,GT = FALSE)

#Identify pK value
FC_pk <- find.pK(sweep.stats = Fem_Coc_Repeated_sweep_stats)
FC_pk <- as.numeric(as.character(FC_pk[which(FC_pk$BCmetric == max(FC_pk$BCmetric)),"pK"]))

#Now model the number of homotypic doublets 
FC_homotypic_props <- modelHomotypic(Fem_Coc_Repeated$CellType)

#Poisson 
nExp_poi <- round(.10*nrow(Fem_Coc_Repeated@meta.data))
nExp_poi_adj <- round(nExp_poi*(1-FC_homotypic_props))

#Find Doublets
#Fem_Coc_Repeated <- doubletFinder_v3(seu = Fem_Coc_Repeated,PCs = 1:17,pN = 0.25,pK = FC_pk,nExp = nExp_poi_adj,reuse.pANN = FALSE,sct = FALSE)
Fem_Coc_Repeated <- doubletFinder_v3(seu = Fem_Coc_Repeated,PCs = 1:17,pN = 0.25,pK = FC_pk,nExp = nExp_poi,reuse.pANN = FALSE,sct = FALSE)
#Change column names
names(Fem_Coc_Repeated@meta.data)[c(16,17)] <- c("pANN","DF.Classification")

#remove all objects except for Fem_Coc_Repeated so that environment does not become overloaded 
rm(FC,FC_metadata,Fem_Coc_Repeated_sweep,Fem_Coc_Repeated_sweep_stats,FC_homotypic_props,FC_pk,nExp_poi,nExp_poi_adj)

##########Male Saline##########
#Pull the Male Saline metadata from the main object that contains the cellnames 
MS_metadata <- subset(All_Groups_Log_Repeated,subset = (Sex_Stim == "Male_Sal"))@meta.data[,c("Non_unique_Cell_Name","CellType","Unique_Cell_Name")]

Male_Sal_Repeated <- All_Groups_log[[3]]

#Male_Sal_Repeated should have rownames that correspond to the non-unique_cell names. MAke that a column and we can merge by that 
Male_Sal_Repeated$Non_unique_Cell_Name <- rownames(Male_Sal_Repeated@meta.data)

#Merge the two dataframes by non-unique-cellname
MS <- merge(x  = Male_Sal_Repeated@meta.data,
            y  = MS_metadata,
            by = "Non_unique_Cell_Name")

#Make non unique names the rownames
rownames(MS) <- MS$Non_unique_Cell_Name

#Make sure nothing is duplicated 
nrow(MS %>% distinct()) == nrow(MS_metadata) #TRUE

#Make sure these are in the same order
#I want MS to be in the same order as Male_Sal_Repeated, therefore MS comes second and Male_Sal_Repeated comes first 
MS <- MS[match(row.names(Male_Sal_Repeated@meta.data),row.names(MS)),]

#Check everything matches
all(row.names(MS) == row.names(Male_Sal_Repeated@meta.data))
all(length(row.names(MS)) == length(row.names(Male_Sal_Repeated@meta.data)))

# #Add the cell type to the metadata
# Male_Sal_Repeated <- AddMetaData(object = Male_Sal_Repeated,metadata = MS$CellType,col.name = "CellType")

#Scale data and dimensionality reduction
Male_Sal_Repeated <- ScaleData(Male_Sal_Repeated,verbose = FALSE)
Male_Sal_Repeated <- RunPCA(Male_Sal_Repeated,npcs = 17,verbose = FALSE)
Male_Sal_Repeated <- RunUMAP(Male_Sal_Repeated,reduction = "pca",dims = 1:17)

#Begin doublet finder workflow with first sweeping parameters 
Male_Sal_Repeated_Sweep <- paramSweep_v3(Male_Sal_Repeated,PCs = 1:17,sct = FALSE)

#summarize the sweep 
Male_Sal_Repeated_Sweep_stats <- summarizeSweep(sweep.list = Male_Sal_Repeated_Sweep,GT = FALSE)

#Run find.pK
MS_pk <- find.pK(sweep.stats = Male_Sal_Repeated_Sweep_stats)

#Identify pK value
MS_pk <- as.numeric(as.character(MS_pk[which(MS_pk$BCmetric == max(MS_pk$BCmetric)),"pK"]))

#Now model the number of homotypic doublets 
MS_homotypic_props <- modelHomotypic(Male_Sal_Repeated$CellType)

#Poisson 
nExp_poi <- round(.10*nrow(Male_Sal_Repeated@meta.data))
nExp_poi_adj <- round(nExp_poi*(1-MS_homotypic_props))

#Find Doublets with adjusted expectation
#Male_Sal_Repeated <- doubletFinder_v3(seu = Male_Sal_Repeated,PCs = 1:17,pN = 0.25,pK = MS_pk,nExp = nExp_poi_adj,reuse.pANN = FALSE,sct = FALSE)
Male_Sal_Repeated <- doubletFinder_v3(seu = Male_Sal_Repeated,PCs = 1:17,pN = 0.25,pK = MS_pk,nExp = nExp_poi,reuse.pANN = FALSE,sct = FALSE)
names(Male_Sal_Repeated@meta.data)[c(16,17)] <- c("pANN","DF.Classification")

#remove all objects except for Male_Sal_Repeated so that environment does not become overloaded 
rm(MS,MS_metadata,MS_homotypic_props,MS_pk,nExp_poi,nExp_poi_adj)

####Male Cocaine#######
#Pull the FeMale_Sal_Repeatedine metadata from the main object that contains the cellnames 
MC_metadata <- subset(All_Groups_Log_Repeated,subset = (Sex_Stim == "Male_Coc"))@meta.data[,c("Non_unique_Cell_Name","CellType","Unique_Cell_Name")]

Male_Coc_Repeated <- All_Groups_log[[4]]

#Fem_Sal_Repeated should have rownames that correspond to the non-unique_cell names. MAke that a column and we can merge by that 
Male_Coc_Repeated$Non_unique_Cell_Name <- rownames(Male_Coc_Repeated@meta.data)

#Merge the two dataframes by non-unique-cellname
MC <- merge(x  = Male_Coc_Repeated@meta.data,
            y  = MC_metadata,
            by = "Non_unique_Cell_Name")

#Make non unique names the rownames
rownames(MC) <- MC$Non_unique_Cell_Name

#Make sure nothing is duplicated 
nrow(MC %>% distinct()) == nrow(MC_metadata) #TRUE

#Now make sure everything is in the same order
MC <- MC[match(row.names(Male_Coc_Repeated@meta.data),row.names(MC)),]

#Check everything matches
all(row.names(MC) == row.names(Male_Coc_Repeated@meta.data))
all(length(row.names(MC)) == length(row.names(Male_Coc_Repeated@meta.data)))

# #Add the cell type to the metadata
# Male_Coc_Repeated <- AddMetaData(object = Male_Coc_Repeated,metadata = MC$CellType,col.name = "CellType")

#Scale data and do dimensionality reduction
Male_Coc_Repeated <- ScaleData(Male_Coc_Repeated,verbose = FALSE)
Male_Coc_Repeated <- RunPCA(Male_Coc_Repeated,npcs = 17,verbose = FALSE)
Male_Coc_Repeated <- RunUMAP(Male_Coc_Repeated,reduction = "pca",dims = 1:17)

#Sweep parameters for optimal pk 
Male_Coc_Repeated_sweep <- paramSweep_v3(Male_Coc_Repeated,PCs = 1:17,sct = FALSE)

#Sweep the above command 
Male_Coc_Repeated_sweep_stats <- summarizeSweep(sweep.list = Male_Coc_Repeated_sweep,GT = FALSE)

#Find the optimal pk value
MC_pk <- find.pK(sweep.stats = Male_Coc_Repeated_sweep_stats)

#Identify pK value
MC_pk <- as.numeric(as.character(MC_pk[which(MC_pk$BCmetric == max(MC_pk$BCmetric)),"pK"]))

#Now model the number of homotypic doublets 
MC_homotypic_props <- modelHomotypic(Male_Coc_Repeated$CellType)

#Poisson 
nExp_poi <- round(.10*nrow(Male_Coc_Repeated@meta.data))
nExp_poi_adj <- round(nExp_poi*(1-MC_homotypic_props))

#Find Doublets with two different stringencies
#Male_Coc_Repeated <- doubletFinder_v3(seu = Male_Coc_Repeated,PCs = 1:17,pN = 0.25,pK = MC_pk,nExp = nExp_poi_adj,reuse.pANN = FALSE,sct = FALSE)
Male_Coc_Repeated <- doubletFinder_v3(seu = Male_Coc_Repeated,PCs = 1:17,pN = 0.25,pK = MC_pk,nExp = nExp_poi,reuse.pANN = FALSE,sct = FALSE)
names(Male_Coc_Repeated@meta.data)[c(16,17)] <- c("pANN","DF.Classification")

rm(MC,MC_metadata,MC_homotypic_props,MC_pk,nExp_poi,nExp_poi_adj)
##############################

#Rbind all of the metadat now. 
Combo_MD <- rbind(Fem_Sal_Repeated@meta.data,Fem_Coc_Repeated@meta.data,Male_Sal_Repeated@meta.data,Male_Coc_Repeated@meta.data)

#Check everything matches
all(Combo_MD$Non_unique_Cell_Name == All_Groups_Log_Repeated$Non_unique_Cell_Name) #TRUE
all(nrow(Combo_MD) == nrow(All_Groups_Log_Repeated@meta.data)) #TRUE

#Add the doublet classification to the metadata with AddMetaData function
All_Groups_Log_Repeated <- AddMetaData(object = All_Groups_Log_Repeated,metadata = Combo_MD$DF.Classification,col.name = "Doublet_Classification")

Idents(All_Groups_Log_Repeated) <- All_Groups_Log_Repeated$Doublet_Classification

Dubs_Umap <- DimPlot(object = All_Groups_Log_Repeated,reduction = "umap")

ggsave(plot = Dubs_Umap,filename = "./MCN_Code/Plots/Repeated_umap_doublets_marked.pdf")

ggsave(plot = cowplot::plot_grid(D_Umap,Dubs_Umap),filename = "./MCN_Code/Plots/Repeated_Decont_Dubs_umaps_sidebyside.pdf")

#Calculate the percent of doublets in each cluster. 
#Create an empty dataframe
Doublet_DF <- data.frame(row.names = levels(All_Groups_Log_Repeated$CellType),
                         Celltype  = levels(All_Groups_Log_Repeated$CellType),
                         Doublet   = NA,
                         Singlet   = NA)
#Add values to the dataframe
for(i in row.names(Doublet_DF)){
  print(i)
  Doublet_DF[i,"Doublet"] <- as.numeric(table(subset(All_Groups_Log_Repeated,subset=(CellType == i))$Doublet_Classification)["Doublet"])
  Doublet_DF[i,"Singlet"] <- as.numeric(table(subset(All_Groups_Log_Repeated,subset=(CellType == i))$Doublet_Classification)["Singlet"])
}
#Calculate total cells
Doublet_DF$Total_Cells <- Doublet_DF$Doublet + Doublet_DF$Singlet
#Now go ahead and calculate percentages
Doublet_DF$Doublet_Pct <- (Doublet_DF$Doublet/Doublet_DF$Total_Cells)*100
Doublet_DF$Singlet_Pct <- (Doublet_DF$Singlet/Doublet_DF$Total_Cells)*100
#Make a bargraph of doublet and singlet percentages
Dub_DF_melt <- reshape2::melt(Doublet_DF)
Dub_DF_melt <- subset(Dub_DF_melt,subset=(variable == "Doublet_Pct" | variable == "Singlet_Pct"))
Dub_DF_melt$variable <- ifelse(Dub_DF_melt$variable == "Doublet_Pct",
                               "Doublet",
                               "Singlet")

Dub_DF_melt$Celltype <- factor(x = Dub_DF_melt$Celltype,
                               levels = unique(Dub_DF_melt$Celltype))

Dub_bargraph <- ggplot(data = Dub_DF_melt,aes(x = Celltype, y = value, fill = variable)) +
  geom_bar(position = "stack",stat = "identity") +
  theme_bw() +
  NoGrid() +
  labs(x = "Cell Type",
       y = "% of Cluster",
       fill = "Cell Class") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "./MCN_Code/Plots/Repeated_Doublets_bargraph.pdf",plot = Dub_bargraph)

#Save the object with doublets identified. 
saveRDS(object = All_Groups_Log_Repeated,
        file   = "./MCN_Code/Objects/Repeated_Cocaine_DubsIdentified.RDS")


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
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ROCR_1.0-11         KernSmooth_2.23-17  fields_11.6         spam_2.5-1          dotCall64_1.0-0    
# [6] dplyr_1.0.2         DoubletFinder_2.0.3 ggplot2_3.3.2       SoupX_1.5.2         SeuratObject_4.0.2 
# [11] Seurat_4.0.4       
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
# [41] reshape2_1.4.4              maps_3.3.0                  rappdirs_0.3.1              Rcpp_1.0.7                 
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
# [85] lattice_0.20-41             purrr_0.3.4                 tensor_1.5                  GenomicAlignments_1.24.0   
# [89] patchwork_1.0.1             htmlwidgets_1.5.2           labeling_0.3                cowplot_1.1.0              
# [93] tidyselect_1.1.0            RcppAnnoy_0.0.19            plyr_1.8.6                  magrittr_1.5               
# [97] R6_2.4.1                    IRanges_2.22.2              generics_0.0.2              DelayedArray_0.14.1        
# [101] pillar_1.4.6                withr_2.3.0                 mgcv_1.8-33                 fitdistrplus_1.1-1         
# [105] survival_3.2-7              abind_1.4-5                 RCurl_1.98-1.2              tibble_3.0.4               
# [109] future.apply_1.6.0          crayon_1.3.4                spatstat.geom_2.4-0         plotly_4.9.2.1             
# [113] data.table_1.13.0           digest_0.6.26               xtable_1.8-4                tidyr_1.1.2                
# [117] httpuv_1.5.4                stats4_4.0.2                munsell_0.5.0               viridisLite_0.3.0
# 
