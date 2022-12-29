#Set working directory and seed
setwd("/data/project/daylab/2019-JD-0040/MCN_Code")
set.seed(1234)

######Load libraries
library(Seurat)
library(babelgene)
library(pheatmap)
library(RColorBrewer)

#Load the NAc_Combo object
NAc_Combo <- readRDS(file = "Objects/NAc_Combo_Integrated.RDS")

DimPlot(object = NAc_Combo,reduction = "umap",label = TRUE) + NoLegend()

#Subset for only neurons
NAc_Neurons <- subset(NAc_Combo,idents = c("Drd1-MSN-1","Drd1-MSN-2","Drd2-MSN-1","Drd2-MSN-2",
                                           "Drd3-MSN","Grm8-MSN","GABAergic","Chat-Interneuron",
                                           "Pvalb-Interneuron","Sst-Interneuron","Glutamatergic"))

#Change the hyphen to a period
Idents(NAc_Neurons) <- gsub(Idents(NAc_Neurons),pattern = "-",replacement = ".")

#Pull log normalized count values
logCounts <- as.data.frame(AverageExpression(object = NAc_Neurons,assays = "RNA",slot = "data")$RNA)

#Make a gene name column
logCounts$GeneName <- rownames(logCounts)

#Make a list consisting of every DEG table from the DESeq2 analysis for Celltype
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

#Now map the human genes to monkey
Monkey_Orthos <- orthologs(genes = DEGs_df$Human_Gene_Symbol,species = "Macaca mulatta")
colnames(Monkey_Orthos)[5] <- "Monkey_Symbol"

#Now merge the DEGs_df and the monkey symbol. 
DEGs_df_all_species <-merge(x = DEGs_df,
                            y = Monkey_Orthos,
                            by.x = "Human_Gene_Symbol",
                            by.y = "human_symbol") #2332 DEGs

#Now load the NHP seurat object
NHP_MSNs <- readRDS(file = "Results_MSNs_processed_final.rds")

#Set identities to MSN_type
Idents(NHP_MSNs) <- NHP_MSNs$MSN_type

#Plot the umap
DimPlot(object = NHP_MSNs,reduction = "umap",label = TRUE) + NoLegend()

###Generate the matrix for rats
#t stat matrix
t_stat_mat <- matrix(nrow = length(unique(DEGs_df_all_species$GeneName)),ncol = length(levels(NAc_Neurons)))

#Change the row names and col names
rownames(t_stat_mat) <- unique(DEGs_df_all_species$GeneName)
colnames(t_stat_mat) <- levels(NAc_Neurons)

for(i in rownames(t_stat_mat)){
  print(i)
  write(i,file="Rn7_Nac_d_t_Calcs_Progress.txt",append=TRUE)
  for(l in colnames(t_stat_mat)){
    #First get all other celltypes
    All_Others <- levels(NAc_Neurons)[which(levels(NAc_Neurons) != l)]
    #Pull gene expression values for cells in cluster of interest
    c1_vals <- GetAssayData(object = NAc_Neurons,assay = "RNA",slot = "data")[i,WhichCells(NAc_Neurons,idents = l)]
    #Pull gene expression values of all other cells
    c2_vals <- GetAssayData(object = NAc_Neurons,assay = "RNA",slot = "data")[i,WhichCells(NAc_Neurons,idents = All_Others)]
    #Get the mean for cell type of interest
    x_c1  <- mean(c1_vals)
    #Get the mean for all other cell types
    x_c2 <- mean(c2_vals)
    #Get standard deviation of c1 vals
    sdev_c1 <- sd(x = c1_vals)
    #Get standard deviation of c1c cals
    sdec_c2 <- sd(x = c2_vals)
    #Now calculate t
    t_stat_mat[i,l] <- (x_c1-x_c2)/sqrt((sdev_c1^2+sdec_c2^2)/2)*sqrt(length(WhichCells(NAc_Neurons,idents = l)))
  }
}

#Write out the t_stat_mat for rats
write.table(x         = t_stat_mat,
            file      = "Rn7_Neuron_DEG_Orthologs_t_stat_mat.txt",
            sep       = "\t",
            col.names = TRUE,
            row.names = TRUE,
            quote     = FALSE)

#Pull the genes that are expressed
NHP_mat_all <- GetAssayData(object = NHP_MSNs,assay = "RNA",slot = "data")
NHP_mat_expressed <- NHP_mat_all[which(rownames(NHP_mat_all) %in% unique(DEGs_df_all_species$Monkey_Symbol)),]

#Now build the same matrix but for NHP
t_stat_mat_NHP <- matrix(nrow = nrow(NHP_mat_expressed),ncol = length(levels(NHP_MSNs)))

#Change the row names and col names
rownames(t_stat_mat_NHP) <- rownames(NHP_mat_expressed)
colnames(t_stat_mat_NHP) <- levels(NHP_MSNs)

for(i in rownames(t_stat_mat_NHP)){
  print(i)
  write(i,file="NHP_MSNs_d_t_Calcs_Progress.txt",append=TRUE)
  for(l in colnames(t_stat_mat_NHP)){
    #First get all other celltypes
    All_Others <- levels(NHP_MSNs)[which(levels(NHP_MSNs) != l)]
    #Pull gene expression values for cells in cluster of interest
    c1_vals <- NHP_mat_expressed[i,WhichCells(NHP_MSNs,idents = l)]
    #Pull gene expression values of all other cells
    c2_vals <- NHP_mat_expressed[i,WhichCells(NHP_MSNs,idents = All_Others)]
    #Get the mean for cell type of interest
    x_c1  <- mean(c1_vals)
    #Get the mean for all other cell types
    x_c2 <- mean(c2_vals)
    #Get standard deviation of c1 vals
    sdev_c1 <- sd(x = c1_vals)
    #Get standard deviation of c1c cals
    sdec_c2 <- sd(x = c2_vals)
    #Now calculate t
    t_stat_mat_NHP[i,l] <- (x_c1-x_c2)/sqrt((sdev_c1^2+sdec_c2^2)/2)*sqrt(length(WhichCells(NHP_MSNs,idents = l)))
  }
}

#Write out the t_stat_mat_NHP for rats
write.table(x         = t_stat_mat_NHP,
            file      = "NHP_MSNs_DEG_Orthologs_t_stat_mat.txt",
            sep       = "\t",
            col.names = TRUE,
            row.names = TRUE,
            quote     = FALSE)


#Now I need to put everything in the same order.
#Screwed up the matrices and now need to load them back in. 
t_stat_mat_rat <- read.table(file = "Rn7_Neuron_DEG_Orthologs_t_stat_mat.txt",sep = "\t",header = TRUE)
colnames(t_stat_mat_rat) <- paste0(colnames(t_stat_mat_rat),"_Rat")
t_stat_mat_rat$gene <- rownames(t_stat_mat_rat)
#now nhp
t_stat_mat_NHP <- read.table(file = "NHP_MSNs_DEG_Orthologs_t_stat_mat.txt",sep = "\t",header = TRUE)
t_stat_mat_NHP <- as.data.frame(t_stat_mat_NHP)
colnames(t_stat_mat_NHP) <- paste0(colnames(t_stat_mat_NHP),"_NHP")
t_stat_mat_NHP$gene <- rownames(t_stat_mat_NHP)

#Merge with gene name key to get the monkey symbol in the dataframe
t_stat_mat_rat <- merge(x = t_stat_mat_rat,
                        y = DEGs_df_all_species[,c("GeneName","Monkey_Symbol")],
                        by.x = "gene",
                        by.y = "GeneName")

#Now merge the two datafrae
t_stat_mat_all <- merge(x    = t_stat_mat_rat,
                        y    = t_stat_mat_NHP,
                        by.x = "Monkey_Symbol",
                        by.y = "gene")


Rat_NHP_Cor_mat <- cor(t_stat_mat_all[,grep(colnames(t_stat_mat_all),pattern = "_Rat")],t_stat_mat_all[,grep(colnames(t_stat_mat_all),pattern = "_NHP")])

colrange <-  seq(-.60,.60, by = 0.01)
colorpal <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(colrange))


#Rearrange the matrix
Rat_NHP_Cor_mat <- Rat_NHP_Cor_mat[c("Drd1.MSN.1_Rat","Drd1.MSN.2_Rat","Drd2.MSN.1_Rat","Drd2.MSN.2_Rat",
                                     "Drd3.MSN_Rat","Grm8.MSN_Rat","Chat.Interneuron_Rat","Pvalb.Interneuron_Rat",
                                     "Sst.Interneuron_Rat","GABAergic_Rat","Glutamatergic_Rat"),
                                   c("D1.Matrix_NHP","D1.Striosome_NHP","D1.Shell.OT_NHP","D1.NUDAP_NHP","D1.ICj_NHP",
                                     "D2.Matrix_NHP","D2.Striosome_NHP","D2.Shell.OT_NHP","D1.D2.Hybrid_NHP")]

pdf(file = "Rat_NHP_Correlations_maxpt6colvalues.pdf")
pheatmap(Rat_NHP_Cor_mat,
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
         legend_breaks=c(seq(-.6,.6, by = 0.3)))
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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] RColorBrewer_1.1-2 pheatmap_1.0.12    babelgene_22.9     SeuratObject_4.0.2 Seurat_4.0.4      
# 
# loaded via a namespace (and not attached):
#   [1] nlme_3.1-149          matrixStats_0.57.0    spatstat.sparse_2.0-0 RcppAnnoy_0.0.19      httr_1.4.2           
# [6] sctransform_0.3.2     tools_4.0.2           R6_2.4.1              irlba_2.3.3           rpart_4.1-15         
# [11] KernSmooth_2.23-17    uwot_0.1.10           mgcv_1.8-33           lazyeval_0.2.2        colorspace_1.4-1     
# [16] tidyselect_1.1.0      gridExtra_2.3         compiler_4.0.2        cli_3.3.0             plotly_4.9.2.1       
# [21] labeling_0.3          scales_1.1.1          lmtest_0.9-38         spatstat.data_2.1-0   ggridges_0.5.2       
# [26] pbapply_1.4-3         rappdirs_0.3.1        goftest_1.2-2         stringr_1.4.0         digest_0.6.26        
# [31] spatstat.utils_2.2-0  pkgconfig_2.0.3       htmltools_0.5.2       fastmap_1.1.0         htmlwidgets_1.5.2    
# [36] rlang_1.0.2           rstudioapi_0.11       shiny_1.5.0           farver_2.0.3          generics_0.0.2       
# [41] zoo_1.8-8             jsonlite_1.7.1        ica_1.0-2             dplyr_1.0.2           magrittr_1.5         
# [46] patchwork_1.0.1       Matrix_1.3-4          Rcpp_1.0.7            munsell_0.5.0         abind_1.4-5          
# [51] reticulate_1.16       lifecycle_0.2.0       stringi_1.7.6         yaml_2.2.1            MASS_7.3-53          
# [56] Rtsne_0.15            plyr_1.8.6            grid_4.0.2            parallel_4.0.2        listenv_0.8.0        
# [61] promises_1.1.1        ggrepel_0.8.2         crayon_1.3.4          miniUI_0.1.1.1        deldir_1.0-6         
# [66] lattice_0.20-41       cowplot_1.1.0         splines_4.0.2         tensor_1.5            pillar_1.4.6         
# [71] igraph_1.2.6          spatstat.geom_2.4-0   future.apply_1.6.0    reshape2_1.4.4        codetools_0.2-16     
# [76] leiden_0.3.3          glue_1.6.2            data.table_1.13.0     png_0.1-7             vctrs_0.4.1          
# [81] httpuv_1.5.4          gtable_0.3.0          RANN_2.6.1            purrr_0.3.4           spatstat.core_2.3-0  
# [86] polyclip_1.10-0       tidyr_1.1.2           scattermore_0.7       future_1.19.1         ggplot2_3.3.2        
# [91] mime_0.9              xtable_1.8-4          later_1.1.0.1         survival_3.2-7        viridisLite_0.3.0    
# [96] tibble_3.0.4          cluster_2.1.0         globals_0.13.1        fitdistrplus_1.1-1    ellipsis_0.3.2       
# [101] ROCR_1.0-11          
# 
# 
# 
