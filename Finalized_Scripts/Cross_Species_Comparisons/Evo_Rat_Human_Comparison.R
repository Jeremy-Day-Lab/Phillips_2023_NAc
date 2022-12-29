#Set working directory and seed
setwd("/data/project/daylab/2019-JD-0040/MCN_Code")
set.seed(1234)

######Load libraries
library(Seurat)
library(babelgene)
library(pheatmap)
library(RColorBrewer)
library(SingleCellExperiment)

#Load the NAc_Combo object
NAc_Combo <- readRDS(file = "Objects/NAc_Combo_Integrated.RDS")

DimPlot(object = NAc_Combo,reduction = "umap",label = TRUE) + NoLegend()

#Subset for only neurons
NAc_Neurons <- subset(NAc_Combo,idents = c("Drd1-MSN-1","Drd1-MSN-2","Drd2-MSN-1","Drd2-MSN-2",
                                           "Drd3-MSN","Grm8-MSN","GABAergic","Chat-Interneuron",
                                           "Pvalb-Interneuron","Sst-Interneuron","Glutamatergic"))

#Change the hyphen to a period
Idents(NAc_Neurons) <- gsub(Idents(NAc_Neurons),pattern = "-",replacement = ".")
Idents(NAc_Combo)   <- gsub(Idents(NAc_Combo),pattern = "-",replacement = ".")

#Now load in the Lieber institute data
load("SCE_NAc-n8_tran-etal.rda")

#Rename the object and remove the old object.
Human_NAc <- sce.nac.tran
rm(sce.nac.tran)

#Pull metadata so you can look at it. 
human_metadata <- colData(Human_NAc)

#Make an object that contains neurons and glia that we want as well as just neurons
#Start with 20571 cells
Human_NAc <- Human_NAc[,Human_NAc$cellType %in% c("Astro_A","Astro_B", #Astrocytes
                                                  "Inhib_A","Inhib_B","Inhib_C","Inhib_D","Inhib_E", #GABAergic neurons? 
                                                  "Macrophage", 
                                                  "Micro","Micro_resting", #Microglia
                                                  "MSN.D1_A","MSN.D1_B","MSN.D1_C","MSN.D1_D","MSN.D1_E","MSN.D1_F", #Drd1-MSNs
                                                  "MSN.D2_A","MSN.D2_B","MSN.D2_C","MSN.D2_D", #Drd2-MSNs
                                                  "Oligo_A","Oligo_B", #oligodendrocytes
                                                  "OPC","OPC_COP")] #Oligo precursor cells
#End with 19892 cells

#Now susbet for just neurons
#Start with 19892 cells
Human_neurons <- Human_NAc[,Human_NAc$cellType %in% c("Inhib_A","Inhib_B","Inhib_C","Inhib_D","Inhib_E",
                                                      "MSN.D1_A","MSN.D1_B","MSN.D1_C","MSN.D1_D","MSN.D1_E","MSN.D1_F", 
                                                      "MSN.D2_A","MSN.D2_B","MSN.D2_C","MSN.D2_D"),] 
#End with 12575 cells

#Read in all of the Rat DEGs. 
#Going to just start with neurons and can always go back adn do entire dataset.
#Make a list consisting of every DEG table from the DESeq2 analysis for Celltype
#Make a gene name column
logCounts <- as.data.frame(AverageExpression(object = NAc_Neurons,assays = "RNA",slot = "data")$RNA)
logCounts$GeneName <- rownames(logCounts)

DEG_Lists <- vector(mode = "list",length = length(levels(Idents(NAc_Neurons))))

#Change names of the list
names(DEG_Lists) <- levels(Idents(NAc_Neurons))

#Read in the DEG lists
for(i in names(DEG_Lists)){
  DEG_Lists[[i]] <- subset(read.table(file = paste0("Tables/DESeq2_Neurons/",i,".txt"),header = TRUE),subset=(padj <= 0.05 & log2FoldChange > 0 & Pct_CellType_Expressing >= 10))
  DEG_Lists[[i]] <- DEG_Lists[[i]][order(DEG_Lists[[i]]$log2FoldChange,decreasing = TRUE)[1:250],]
  DEG_Lists[[i]] <- merge(x  = DEG_Lists[[i]],
                          y  = logCounts[,c(i,"GeneName")],
                          by = "GeneName")
  colnames(DEG_Lists[[i]])[10] <- "logCount"
  DEG_Lists[[i]]$CellType <- i
}

#Create a dataframe
DEGs_df <- do.call(what = rbind, DEG_Lists)

#Remove some unneeded columns
DEGs_df <- DEGs_df[,-which(colnames(DEGs_df) %in% c("Pct_Expressing","Pct_CellType_Expressing","Pct_Other_Expressing"))]

#Read in the homology report from JAX lab
hom <- read.delim(file   = "MAGMA/Jax_homology_Report_AllOrganisms_111522.txt",
                  header = TRUE,
                  sep    = "\t")

#Pull human and rat
hom_hs  <- hom[hom$Common.Organism.Name == "human", ]
hom_rat <- hom[hom$Common.Organism.Name == "rat", ]

#Create new columns for the dataframes 
DEGs_df$Human_Gene_Symbol <- NA
#Now map to the human homolog
#for the whole gene set
for(i in 1:nrow(DEGs_df)){
  print(i)
  #Only run the code if a homolog is found. 
  if(length(hom_hs[hom_hs$DB.Class.Key %in% hom_rat[hom_rat$Symbol %in% DEGs_df[i,"GeneName"],"DB.Class.Key"],"Symbol"]) > 0){
    DEGs_df[i,"Human_Gene_Symbol"] <- hom_hs[hom_hs$DB.Class.Key %in% hom_rat[hom_rat$Symbol %in% DEGs_df[i,"GeneName"],"DB.Class.Key"],"Symbol"][1]
  }else{
    next
  }
}

#Remove NA values
DEGs_df <- na.omit(DEGs_df) #start with 2750. End with 2379. Not terrible. 

#Pull the genes that are expressed
human_mat_all <- assay(Human_neurons,"logcounts")
human_mat_expressed <- human_mat_all[which(rownames(human_mat_all) %in% unique(DEGs_df$Human_Gene_Symbol)),] 

#Now build the same matrix but for human
t_stat_mat_human <- matrix(nrow = nrow(human_mat_expressed),ncol = length(unique(Human_neurons$cellType)))

#Change the row names and col names
rownames(t_stat_mat_human) <- rownames(human_mat_expressed)
colnames(t_stat_mat_human) <- unique(Human_neurons$cellType)


for(i in rownames(t_stat_mat_human)){
  print(i)
  for(l in colnames(t_stat_mat_human)){
    #First get all other celltypes
    All_Others <- colnames(t_stat_mat_human)[which(colnames(t_stat_mat_human) != l)]
    #Pull gene expression values for cells in cluster of interest
    c1_vals <- human_mat_expressed[i,rownames(colData(Human_neurons)[which(Human_neurons$cellType %in% colnames(t_stat_mat_human)[which(colnames(t_stat_mat_human) == l)]),])]
    #Pull gene expression values of all other cells
    c2_vals <-  human_mat_expressed[i,rownames(colData(Human_neurons)[which(Human_neurons$cellType %in% colnames(t_stat_mat_human)[which(colnames(t_stat_mat_human) != l)]),])]
    #Get the mean for cell type of interest
    x_c1  <- mean(c1_vals)
    #Get the mean for all other cell types
    x_c2 <- mean(c2_vals)
    #Get standard deviation of c1 vals
    sdev_c1 <- sd(x = c1_vals)
    #Get standard deviation of c1c cals
    sdec_c2 <- sd(x = c2_vals)
    #Now calculate t
    t_stat_mat_human[i,l] <- (x_c1-x_c2)/sqrt((sdev_c1^2+sdec_c2^2)/2)*sqrt(length(rownames(colData(Human_neurons)[which(Human_neurons$cellType %in% colnames(t_stat_mat_human)[which(colnames(t_stat_mat_human) == l)]),])))
  }
}

#Write out the t_stat_mat_NHP for rats
write.table(x         = t_stat_mat_human,
            file      = "human_neurons_DEG_Orthologs_t_stat_mat.txt",
            sep       = "\t",
            col.names = TRUE,
            row.names = TRUE,
            quote     = FALSE)


#Now read in the data. 
t_stat_mat_rat <- read.table(file = "Rn7_Neuron_DEG_Orthologs_t_stat_mat.txt",sep = "\t")
colnames(t_stat_mat_rat) <- paste0(colnames(t_stat_mat_rat),"_Rat")
t_stat_mat_rat$gene <- rownames(t_stat_mat_rat)

t_stat_mat_rat <- merge(x = t_stat_mat_rat,
                        y = DEGs_df[,c("GeneName","Human_Gene_Symbol")],
                        by.x = "gene",
                        by.y = "GeneName")

#Change colnames for human
t_stat_mat_human <- as.data.frame(t_stat_mat_human)
colnames(t_stat_mat_human) <- paste0(colnames(t_stat_mat_human),"_human")
t_stat_mat_human$gene <- rownames(t_stat_mat_human)

#Some values in the t_stat_human_matrix are NAs because the expression values are 0s. Need to remove that. 
t_stat_mat_human <- t_stat_mat_human[-which(is.na(t_stat_mat_human)),]

#Now merge the two datafrae
t_stat_mat_all <- merge(x    = t_stat_mat_rat,
                        y    = t_stat_mat_human,
                        by.x = "Human_Gene_Symbol",
                        by.y = "gene") #2282

Rat_Human_Cor_mat <- cor(t_stat_mat_all[,grep(colnames(t_stat_mat_all),pattern = "_Rat")],
                         t_stat_mat_all[,grep(colnames(t_stat_mat_all),pattern = "_human")])

colrange <-  seq(-.65,.65, by = 0.01)
colorpal <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(colrange))


#Rearrange the matrix
Rat_Human_Cor_mat <- Rat_Human_Cor_mat[c("Drd1.MSN.1_Rat","Drd1.MSN.2_Rat","Drd2.MSN.1_Rat","Drd2.MSN.2_Rat",
                                         "Drd3.MSN_Rat","Grm8.MSN_Rat","Chat.Interneuron_Rat","Pvalb.Interneuron_Rat",
                                         "Sst.Interneuron_Rat","GABAergic_Rat","Glutamatergic_Rat"),
                                       colnames(Rat_Human_Cor_mat)[order(colnames(Rat_Human_Cor_mat))]]

pdf(file = "Rat_Human_Neurons_Correlations.pdf")
pheatmap(Rat_Human_Cor_mat,
         color=colorpal,
         cluster_cols=F, 
         cluster_rows=F,
         breaks=colrange,
         fontsize=11, 
         fontsize_row=11.5, 
         fontsize_col=12,
         display_numbers=T, 
         number_format="%.2f", 
         fontsize_number=6.5,
         legend_breaks=c(seq(-.65,.65, by = 0.325)))
dev.off()

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
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] SingleCellExperiment_1.10.1 SummarizedExperiment_1.18.2 DelayedArray_0.14.1         matrixStats_0.57.0         
# [5] Biobase_2.48.0              GenomicRanges_1.40.0        GenomeInfoDb_1.24.2         IRanges_2.22.2             
# [9] S4Vectors_0.26.1            BiocGenerics_0.34.0         RColorBrewer_1.1-2          pheatmap_1.0.12            
# [13] babelgene_22.9              SeuratObject_4.0.2          Seurat_4.0.4               
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15             colorspace_1.4-1       deldir_1.0-6           ellipsis_0.3.2         ggridges_0.5.2        
# [6] XVector_0.28.0         rstudioapi_0.11        spatstat.data_2.1-0    farver_2.0.3           leiden_0.3.3          
# [11] listenv_0.8.0          ggrepel_0.8.2          codetools_0.2-16       splines_4.0.2          polyclip_1.10-0       
# [16] jsonlite_1.7.1         ica_1.0-2              cluster_2.1.0          png_0.1-7              uwot_0.1.10           
# [21] shiny_1.5.0            sctransform_0.3.2      spatstat.sparse_2.0-0  compiler_4.0.2         httr_1.4.2            
# [26] Matrix_1.3-4           fastmap_1.1.0          lazyeval_0.2.2         cli_3.3.0              later_1.1.0.1         
# [31] htmltools_0.5.2        tools_4.0.2            igraph_1.2.6           gtable_0.3.0           glue_1.6.2            
# [36] GenomeInfoDbData_1.2.3 RANN_2.6.1             reshape2_1.4.4         dplyr_1.0.2            rappdirs_0.3.1        
# [41] Rcpp_1.0.7             scattermore_0.7        vctrs_0.4.1            nlme_3.1-149           lmtest_0.9-38         
# [46] stringr_1.4.0          globals_0.13.1         mime_0.9               miniUI_0.1.1.1         lifecycle_0.2.0       
# [51] irlba_2.3.3            goftest_1.2-2          future_1.19.1          MASS_7.3-53            zlibbioc_1.34.0       
# [56] zoo_1.8-8              scales_1.1.1           spatstat.core_2.3-0    promises_1.1.1         spatstat.utils_2.2-0  
# [61] yaml_2.2.1             reticulate_1.16        pbapply_1.4-3          gridExtra_2.3          ggplot2_3.3.2         
# [66] rpart_4.1-15           stringi_1.7.6          rlang_1.0.2            pkgconfig_2.0.3        bitops_1.0-6          
# [71] lattice_0.20-41        ROCR_1.0-11            purrr_0.3.4            tensor_1.5             labeling_0.3          
# [76] patchwork_1.0.1        htmlwidgets_1.5.2      cowplot_1.1.0          tidyselect_1.1.0       RcppAnnoy_0.0.19      
# [81] plyr_1.8.6             magrittr_1.5           R6_2.4.1               generics_0.0.2         pillar_1.4.6          
# [86] mgcv_1.8-33            fitdistrplus_1.1-1     survival_3.2-7         abind_1.4-5            RCurl_1.98-1.2        
# [91] tibble_3.0.4           future.apply_1.6.0     crayon_1.3.4           KernSmooth_2.23-17     spatstat.geom_2.4-0   
# [96] plotly_4.9.2.1         grid_4.0.2             data.table_1.13.0      digest_0.6.26          xtable_1.8-4          
# [101] tidyr_1.1.2            httpuv_1.5.4           munsell_0.5.0          viridisLite_0.3.0     
# 
# 
# 
# 
# 
# 
# 
