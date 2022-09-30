#Calculate cell type specific genes. To do this you need to create a single matrix where the samples are the Celltype. 
#Set working directory and seed
setwd("/data/project/daylab/")
set.seed(1234)

######Load libraries
library(Seurat)
library(ggplot2)
library(Libra)
library(dplyr)
library(ggplot2)
library(stringr)
library(ComplexUpset)
library(Matrix.utils)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)

#Load the NAc_Combo object
NAc_Combo <- readRDS(file = "2019-JD-0040/MCN_Code/Objects/NAc_Combo_Integrated.RDS")

DimPlot(object = NAc_Combo,reduction = "umap",label = TRUE) + NoLegend()

#Create a sample id column that is Dataset_Sex_Stim
NAc_Combo$sample.id <- as.factor(paste(NAc_Combo$Dataset,NAc_Combo$Sex_Stim,sep = "_"))

#Get cell and sample metrics for aggregation
NAc_Combo@meta.data$CellType <- as.factor(NAc_Combo@meta.data$Combo_CellType)
cell_names <- purrr::set_names(levels(NAc_Combo@meta.data$CellType))
cell_names
# Polydendrocyte              Olig-1               Mural           Microglia           Astrocyte       Glutamatergic     Sst-Interneuron   Pvalb-Interneuron    Chat-Interneuron           GABAergic 
# "Polydendrocyte"            "Olig-1"             "Mural"         "Microglia"         "Astrocyte"     "Glutamatergic"   "Sst-Interneuron" "Pvalb-Interneuron"  "Chat-Interneuron"         "GABAergic" 
# Grm8-MSN            Drd3-MSN          Drd2-MSN-2          Drd2-MSN-1          Drd1-MSN-2          Drd1-MSN-1 
# "Grm8-MSN"          "Drd3-MSN"        "Drd2-MSN-2"        "Drd2-MSN-1"        "Drd1-MSN-2"        "Drd1-MSN-1" 

#Number of clusters
cluster_total <- length(cell_names)
cluster_total
# [1] 16

# Named vector of sample names
NAc_Combo$Sex_Stim <- factor(NAc_Combo$Sex_Stim)
sample_ids <- purrr::set_names(levels(NAc_Combo@meta.data$Sex_Stim))
sample_ids
# Fem_Coc    Fem_Sal   Male_Coc   Male_Sal 
# "Fem_Coc"  "Fem_Sal" "Male_Coc" "Male_Sal" 

# Total number of samples 
sample_total <- length(sample_ids)
sample_total
# [1] 4

#Figure out how many cells in each main group
table(NAc_Combo@meta.data$Sex_Stim)
# Fem_Coc  Fem_Sal Male_Coc Male_Sal 
# 11432     9407     8921     9494 

##########Count aggregation to sample level
groups <- NAc_Combo@meta.data[, c("CellType", "sample.id")]

# Aggregate across cluster-sample groups (raw counts, thus counts slot)
count_aggr <- Matrix.utils::aggregate.Matrix(t(NAc_Combo@assays$RNA@counts), 
                                             groupings = groups, fun = "sum") 

dim(count_aggr)
# [1]   128 30560


#Transpose count_aggr
count_aggr_t <- t(count_aggr)

#change the count_aggr_t
colnames(count_aggr_t) <- gsub(x = colnames(count_aggr_t),pattern = "-",replacement = ".")

#Create a metadata dataframe
metadata <- data.frame(cluster_id = sub("_.*","", colnames(count_aggr_t)),
                       sample_id  = colnames(count_aggr_t),
                       Sex.Stim  = sub("^[^_]*_","", colnames(count_aggr_t)))
metadata$dataset <- as.factor(as.character(lapply(strsplit(metadata$Sex.Stim,"_"),"[",1)))
metadata$Sex     <- as.factor(as.character(lapply(strsplit(metadata$Sex.Stim,"_"),"[",2)))
metadata$Stim    <- as.factor(as.character(lapply(strsplit(metadata$Sex.Stim,"_"),"[",3)))
rownames(metadata) <- metadata$sample_id

###Calculate DEGs 
counts_mat <- NAc_Combo@assays$RNA@counts
Idents(NAc_Combo) <- gsub(x = Idents(NAc_Combo),pattern = "-",replacement = ".")
for(i in unique(metadata$cluster_id)){
  ############Make the DESEq object############
  metadata$CellType <- ifelse(metadata$cluster_id == i,
                              i,
                              "Other")
  dds <- DESeqDataSetFromMatrix(count_aggr_t, 
                                colData = metadata, 
                                design = ~ dataset + Stim + CellType)
  dds$CellType <- relevel(dds$CellType, ref = "Other")
  keep <- rowMeans(counts(dds)) > 5
  dds <- dds[keep,]
  print(paste("dds made for",i))
  ############Quality control############
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  #PCA for major celltype
  data <- plotPCA(vsd, intgroup = c("CellType"), returnData = TRUE)
  data$cluster_id <- as.character(lapply(strsplit(x = data$name,split = "_"),"[",1))
  percentVar <- round(100 * attr(data, "percentVar"))
  pca.plot.CellType <- ggplot(data, aes(PC1, PC2, color = CellType)) +
    geom_point(size = 7) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme(text = element_text(size = 20)) +
    theme_bw(base_size = 16)
  #PCA for cluster_id
  pca.plot.cluster_id <- ggplot(data, aes(PC1, PC2, color = cluster_id)) +
    geom_point(size = 7) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme(text = element_text(size = 20)) +
    theme_bw(base_size = 16)
  ggsave(filename = paste0("2019-JD-0040/MCN_Code/Plots/DESeq2_CellType/PCA",i,".pdf"),
         plot     = cowplot::plot_grid(plotlist = list(pca.plot.CellType,pca.plot.cluster_id),ncol = 2),
         width = 12,
         height = 7)
  print(paste("PCA plots made for",i))
  ############Get results############
  dds <- DESeq(dds, test="LRT", reduced = ~ dataset + Stim)
  res <- as.data.frame(results(dds))
  res$GeneName <- rownames(res)
  cells.1    <- WhichCells(NAc_Combo,idents = i)
  cells.2    <- setdiff(x = row.names(NAc_Combo@meta.data),y=cells.1)
  res$Pct_Expressing <- NA
  res$Pct_CellType_Expressing <- (rowSums(x = counts_mat[res$GeneName, cells.1,drop = FALSE] > 0) / length(cells.1))*100
  res$Pct_Other_Expressing <- (rowSums(x = counts_mat[res$GeneName, cells.2,drop = FALSE] > 0) / length(cells.2))*100
  write.table(file      = paste0("2019-JD-0040/MCN_Code/Tables/DESeq2_CellTypes/",i,".txt"),
              x         = res,
              sep       = "\t",
              col.names = TRUE,
              row.names = FALSE,
              quote     = FALSE)
  print(paste("DEG lists written for",i))
  rm(dds,keep,vsd,data,percentVar,pca.plot.CellType,pca.plot.cluster_id,res)
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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] pheatmap_1.0.12             RColorBrewer_1.1-2          DESeq2_1.28.1               SummarizedExperiment_1.18.2 DelayedArray_0.14.1        
# [6] matrixStats_0.57.0          Biobase_2.48.0              GenomicRanges_1.40.0        GenomeInfoDb_1.24.2         IRanges_2.22.2             
# [11] S4Vectors_0.26.1            BiocGenerics_0.34.0         Matrix.utils_0.9.8          Matrix_1.3-4                ComplexUpset_1.3.3         
# [16] stringr_1.4.0               dplyr_1.0.2                 Libra_1.0.0                 ggplot2_3.3.2               SeuratObject_4.0.2         
# [21] Seurat_4.0.4               
# 
# loaded via a namespace (and not attached):
#   [1] blme_1.0-5             plyr_1.8.6             igraph_1.2.6           lazyeval_0.2.2         TMB_1.8.1              splines_4.0.2         
# [7] BiocParallel_1.22.0    listenv_0.8.0          scattermore_0.7        digest_0.6.26          htmltools_0.5.2        lmerTest_3.1-3        
# [13] magrittr_1.5           memoise_1.1.0          tensor_1.5             cluster_2.1.0          ROCR_1.0-11            limma_3.44.3          
# [19] globals_0.13.1         annotate_1.66.0        tester_0.1.7           spatstat.sparse_2.0-0  colorspace_1.4-1       blob_1.2.1            
# [25] rappdirs_0.3.1         ggrepel_0.8.2          crayon_1.3.4           RCurl_1.98-1.2         jsonlite_1.7.1         genefilter_1.70.0     
# [31] lme4_1.1-26            spatstat.data_2.1-0    survival_3.2-7         zoo_1.8-8              glue_1.6.2             polyclip_1.10-0       
# [37] gtable_0.3.0           zlibbioc_1.34.0        XVector_0.28.0         leiden_0.3.3           future.apply_1.6.0     abind_1.4-5           
# [43] scales_1.1.1           DBI_1.1.0              edgeR_3.30.3           miniUI_0.1.1.1         Rcpp_1.0.7             viridisLite_0.3.0     
# [49] xtable_1.8-4           reticulate_1.16        spatstat.core_2.3-0    bit_4.0.4              htmlwidgets_1.5.2      httr_1.4.2            
# [55] ellipsis_0.3.2         ica_1.0-2              farver_2.0.3           pkgconfig_2.0.3        XML_3.99-0.5           uwot_0.1.10           
# [61] deldir_1.0-6           locfit_1.5-9.4         labeling_0.3           tidyselect_1.1.0       rlang_1.0.2            reshape2_1.4.4        
# [67] later_1.1.0.1          AnnotationDbi_1.50.3   pbmcapply_1.5.0        munsell_0.5.0          tools_4.0.2            cli_3.3.0             
# [73] generics_0.0.2         RSQLite_2.2.1          ggridges_0.5.2         fastmap_1.1.0          yaml_2.2.1             goftest_1.2-2         
# [79] bit64_4.0.5            fitdistrplus_1.1-1     purrr_0.3.4            RANN_2.6.1             pbapply_1.4-3          future_1.19.1         
# [85] nlme_3.1-149           mime_0.9               grr_0.9.5              compiler_4.0.2         rstudioapi_0.11        plotly_4.9.2.1        
# [91] png_0.1-7              spatstat.utils_2.2-0   tibble_3.0.4           statmod_1.4.35         geneplotter_1.66.0     stringi_1.7.6         
# [97] forcats_0.5.0          lattice_0.20-41        nloptr_1.2.2.2         vctrs_0.4.1            pillar_1.4.6           lifecycle_0.2.0       
# [103] spatstat.geom_2.4-0    lmtest_0.9-38          RcppAnnoy_0.0.19       data.table_1.13.0      cowplot_1.1.0          bitops_1.0-6          
# [109] irlba_2.3.3            httpuv_1.5.4           patchwork_1.0.1        R6_2.4.1               promises_1.1.1         KernSmooth_2.23-17    
# [115] gridExtra_2.3          codetools_0.2-16       boot_1.3-25            MASS_7.3-53            withr_2.3.0            sctransform_0.3.2     
# [121] GenomeInfoDbData_1.2.3 mgcv_1.8-33            grid_4.0.2             rpart_4.1-15           glmmTMB_1.1.3          tidyr_1.1.2           
# [127] minqa_1.2.4            Rtsne_0.15             numDeriv_2016.8-1.1    shiny_1.5.0   