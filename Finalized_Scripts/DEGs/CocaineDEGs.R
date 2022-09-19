#Cocaine DEGs
setwd("/data/project/daylab/")
set.seed(1234)
#This workflow is based off work from Lara Ianov, Ph.D. Thanks Lara! 
library(Libra)
library(dplyr)
library(Seurat)
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
####Function#####
deter_direc <- function(x) {
  
  log_names <- c("avg_logFC", "avg_log2FC", "log2FC", "log2FoldChange")
  
  log_select <- x %>% select(any_of(log_names))
  
  ifelse(log_select <= 0, "down", ifelse(log_select >= 0, "up", "no_change"))
}

fromList <- function (input) {
  # Same as original UpSetR::fromList(), but modified as shown in https://github.com/hms-dbmi/UpSetR/issues/85
  # thanks to @docmanny
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names! (in this case gene names)
  row.names(data) <- elements
  return(data)
}

plot_output <- function(p, file_name, w_png=700, h_png=600, w_pdf=12, h_pdf=8, show_plot = TRUE){
  
  png(paste0(file_name,".png"), width = w_png, height = h_png)
  plot(eval(p))
  dev.off()
  
  pdf(paste0(file_name,".pdf"), width = w_pdf, height = h_pdf)
  plot(eval(p))
  dev.off()
  
  if (show_plot) {
    plot(eval(plot))
  }
  
}

# for pheatmaps as it uses grid system instead of ggplot
pheatmap_output <- function(x, file_name, w_png=900, h_png=700, w_pdf=12, h_pdf=8) {
  
  png(paste0(file_name,".png"), width = w_png, height = h_png)
  
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
  
  pdf(paste0(file_name,".pdf"), width = w_pdf, height = h_pdf)
  
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

nested_lapply <- function(data, FUN) {
  lapply(data, function(sublist) { lapply(sublist, FUN) })
}

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
#[1] 4

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
#128 individual samples from 16 clusters and 8 individual GEM wells

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
raw_counts_list <- split.data.frame(
  count_aggr,
  factor(gsub("-", ".",cell_names)) # chaning list names from - to . is optional (gsub)
) %>%
  lapply(function(x) {
    magrittr::set_colnames(t(x), gsub("-", ".",rownames(x)))
  })


####Sample level metadata
get_sample_ids <- function(x){
  raw_counts_list[[x]] %>%
    colnames()
}

de_samples <- purrr::map(1:length(cell_names), get_sample_ids) %>%
  unlist()


# Create a data frame with the sample IDs, cluster IDs and condition
metadata <- data.frame(cluster_id = sub("_.*","", de_samples),
                       sample_id = de_samples,
                       Sex.Stim = sub("^[^_]*_","", de_samples))
#Create extra columns in the metadata
metadata$dataset <- as.factor(as.character(lapply(strsplit(metadata$Sex.Stim,"_"),"[",1)))
metadata$Sex     <- as.factor(as.character(lapply(strsplit(metadata$Sex.Stim,"_"),"[",2)))
metadata$Stim    <- as.factor(as.character(lapply(strsplit(metadata$Sex.Stim,"_"),"[",3)))

#Create a DESEqDataSet from the object
dds_all_cells <- mapply(FUN = function(x,y,z) {
  
  #subset metadata by cell type being tested
  
  cluster_metadata <- subset(metadata, cluster_id == z)
  
  # Assign the rownames of the metadata to be the sample IDs
  rownames(cluster_metadata) <- cluster_metadata$sample_id
  
  print(head(cluster_metadata))
  
  # subset the raw counts to cell type being tested
  counts <- raw_counts_list[[y]]
  
  cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
  
  print(head(cluster_counts))
  
  # check that all of the row names of the metadata are the same and in the same order as the cols of the counts
  stopifnot(all(rownames(cluster_metadata) == colnames(cluster_counts)))
  
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ dataset + Stim)
  
  # set reference to Saline
  
  dds$Stim <- relevel(dds$Stim, ref = "Sal")
  
  return(dds)
  
}, x = raw_counts_list, y = 1:length(raw_counts_list), z = gsub("-", ".",names(raw_counts_list)))

###Remove genes with <5 counts
dds_all_cells <- mapply(FUN = function(x, z) {
  
  #----- Counts Pre-filtering based on rowMeans -------
  # here rowMeans is used, but other approaches that can be used are rowSums or min. per sample
  
  message(paste0("Number of genes before pre-filtering ", z, ": ",  nrow(counts(x))))
  
  keep <- rowMeans(counts(x)) > 5
  
  x <- x[keep,]
  
  message(paste0("Number of genes after pre-filtering ", z, ": ",  nrow(counts(x))))
  
  return(x)
  
}, x = dds_all_cells, z = gsub("-", ".",names(dds_all_cells)))

#Check to see how many rows with counts > 5
mapply(FUN = function(x){
  message(paste0(nrow(counts(x))))
},x = dds_all_cells)


#############QC######################
mapply(FUN = function(x, z) {
  vsd <- varianceStabilizingTransformation(x, blind = FALSE)
  vsd
  
  # PCA
  #-----------------Stim---------------------
  data <- plotPCA(vsd, intgroup = c("Stim"), returnData = TRUE)
  percentVar <- round(100 * attr(data, "percentVar"))
  
  pca.plot.stim <- ggplot(data, aes(PC1, PC2, color = Stim)) +
    geom_point(size = 7) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme(text = element_text(size = 20)) +
    theme_bw(base_size = 16)
  
  #-----------------Sex---------------------
  data <- plotPCA(vsd, intgroup = c("Sex"), returnData = TRUE)
  
  pca.plot.sex <- ggplot(data, aes(PC1, PC2, color = Sex)) +
    geom_point(size = 7) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme(text = element_text(size = 20)) +
    theme_bw(base_size = 16)
  
  #-----------------Sex.Stim---------------------
  data <- plotPCA(vsd, intgroup = c("Sex.Stim"), returnData = TRUE)
  
  pca.plot.sex.stim <- ggplot(data, aes(PC1, PC2, color = Sex.Stim)) +
    geom_point(size = 7) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme(text = element_text(size = 20)) +
    theme_bw(base_size = 16)
  
  all_pca <- pca.plot.stim + pca.plot.sex + pca.plot.sex.stim
  
  file_name_pca <- paste0("2019-JD-0040/MCN_Code/Plots/DESeq2_QC/PCA_", z)
  
  plot_output(all_pca, file_name = file_name_pca, 
              w_png = 1800, w_pdf = 20, h_png = 400, h_pdf = 6, show_plot = FALSE)
  
  
  # sample-to-samples distances heatmap:
  # dist computes distance with Euclidean method
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  heatmap_dist <- pheatmap(sampleDistMatrix,
                           clustering_distance_rows = sampleDists,
                           clustering_distance_cols = sampleDists,
                           col = colors
  )
  
  file_name_heatmap <- paste0("2019-JD-0040/MCN_Code/Plots/DESeq2_QC/sample_to_sample_heatmap_", z)
  
  # save output
  pheatmap_output(heatmap_dist, file_name = file_name_heatmap)
  
  
}, x = dds_all_cells, z = gsub("-", ".", names(dds_all_cells)))


###Calculate DEGs
dds_all_cells <- mapply(FUN = function(x, z) {
  
  x <- DESeq(x, test="LRT", reduced = ~ dataset)
  
  # checking coeff by resultsNames(dds):
  
  resultsNames(x)
  
  # Dispersion plot
  
  file_name_dist <- paste0("2019-JD-0040/MCN_Code/Plots/DESeq2_DEGResults/DispPlot_", z, ".pdf")
  
  pdf(file_name_dist)
  plotDispEsts(x)
  dev.off()
  
  return(x)
  
}, x = dds_all_cells, z = gsub("-", ".", names(dds_all_cells)))


#Calculate results
res_all_cells <- mapply(FUN = function(x) {
  
  res <- results(x) # ref level established earlier by relevel
  
  return(res)
  
}, x = dds_all_cells)


#Now calculate percent cells expressing
Idents(NAc_Combo) <- gsub(Idents(NAc_Combo),pattern = "-",replacement = ".")
SeuratCounts <- t(GetAssayData(object = NAc_Combo,slot = "data",assay = "RNA"))
for(i in names(res_all_cells)){
  print(paste("Beginning:",nrow(res_all_cells[[i]])))
  res_all_cells[[i]] <- as.data.frame(res_all_cells[[i]])
  #Create a gene name column
  res_all_cells[[i]]$Gene <- row.names(res_all_cells[[i]])
  #Pull the genes and calculate percent expressing
  Genes_PctExp <- as.data.frame(colMeans(SeuratCounts[WhichCells(object = NAc_Combo,idents = i),res_all_cells[[i]]$Gene]>0)*100)
  #Change the column names
  Genes_PctExp$Gene <- row.names(Genes_PctExp)
  colnames(Genes_PctExp)[1] <- "Pct_Expressing"
  #Add the percent expressing to the dataframe
  res_all_cells[[i]] <- merge(x = res_all_cells[[i]],
                              y = Genes_PctExp,
                              by = "Gene")
  print(paste("End:",nrow(res_all_cells[[i]])))
}


#Write out the files
mapply(FUN = function(x, z) {
  
  file_name <- paste0("2019-JD-0040/MCN_Code/Tables/DESeq2_Coc_Collapsed/DESeq2_res_", z, ".csv")
  
  write.csv(x, file=file_name, row.names = TRUE,quote = FALSE)
  
}, x = res_all_cells, z = gsub("-", ".", names(res_all_cells)))

#Save counts
mapply(FUN = function(x, z) {
  
  file_name_norm <- paste0("2019-JD-0040/MCN_Code/Tables/Normalized_Counts/DESeq2_Coc_Collapsed/normalized_counts_", z, ".csv")
  
  normalized_counts <- as.data.frame(counts(x, normalized = TRUE))
  
  write.csv(normalized_counts, file=file_name_norm, row.names = TRUE,quote = FALSE)
  
}, x = dds_all_cells, z = gsub("-", ".", names(dds_all_cells)))


#Make a table of the number of DEGs
DEGs_DF <- as.data.frame(matrix(ncol = 2,nrow=16))
colnames(DEGs_DF) <- c("CellType","No_DEGs")
rownames(DEGs_DF) <- names(dds_all_cells)
DEGs_DF$CellType  <- names(dds_all_cells)

for(i in names(dds_all_cells)){
  x <- as.data.frame(res_all_cells[[i]])
  DEGs_DF[i,"No_DEGs"] <- nrow(subset(x,subset=(padj <= 0.05)))
}

write.table(x         = DEGs_DF,
            file      = "2019-JD-0040/MCN_Code/Tables/DESeq2_Coc_Collapsed/No_DEGs_Per_CellType.txt",
            sep       = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE)

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
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] pheatmap_1.0.12             RColorBrewer_1.1-2          DESeq2_1.28.1               SummarizedExperiment_1.18.2 DelayedArray_0.14.1         matrixStats_0.57.0         
# [7] Biobase_2.48.0              GenomicRanges_1.40.0        GenomeInfoDb_1.24.2         IRanges_2.22.2              S4Vectors_0.26.1            BiocGenerics_0.34.0        
# [13] Matrix.utils_0.9.8          Matrix_1.3-4                ComplexUpset_1.3.3          stringr_1.4.0               ggplot2_3.3.2               SeuratObject_4.0.2         
# [19] Seurat_4.0.4                dplyr_1.0.2                 Libra_1.0.0                
# 
# loaded via a namespace (and not attached):
#   [1] blme_1.0-5             plyr_1.8.6             igraph_1.2.6           lazyeval_0.2.2         TMB_1.8.1              splines_4.0.2          BiocParallel_1.22.0    listenv_0.8.0         
# [9] scattermore_0.7        digest_0.6.26          htmltools_0.5.2        lmerTest_3.1-3         magrittr_1.5           memoise_1.1.0          tensor_1.5             cluster_2.1.0         
# [17] ROCR_1.0-11            limma_3.44.3           globals_0.13.1         annotate_1.66.0        tester_0.1.7           spatstat.sparse_2.0-0  colorspace_1.4-1       blob_1.2.1            
# [25] rappdirs_0.3.1         ggrepel_0.8.2          crayon_1.3.4           RCurl_1.98-1.2         jsonlite_1.7.1         genefilter_1.70.0      lme4_1.1-26            spatstat.data_2.1-0   
# [33] survival_3.2-7         zoo_1.8-8              glue_1.6.2             polyclip_1.10-0        gtable_0.3.0           zlibbioc_1.34.0        XVector_0.28.0         leiden_0.3.3          
# [41] future.apply_1.6.0     abind_1.4-5            scales_1.1.1           DBI_1.1.0              edgeR_3.30.3           miniUI_0.1.1.1         Rcpp_1.0.7             viridisLite_0.3.0     
# [49] xtable_1.8-4           reticulate_1.16        spatstat.core_2.3-0    bit_4.0.4              htmlwidgets_1.5.2      httr_1.4.2             ellipsis_0.3.2         ica_1.0-2             
# [57] farver_2.0.3           pkgconfig_2.0.3        XML_3.99-0.5           uwot_0.1.10            deldir_1.0-6           locfit_1.5-9.4         labeling_0.3           tidyselect_1.1.0      
# [65] rlang_1.0.2            reshape2_1.4.4         later_1.1.0.1          AnnotationDbi_1.50.3   pbmcapply_1.5.0        munsell_0.5.0          tools_4.0.2            cli_3.3.0             
# [73] generics_0.0.2         RSQLite_2.2.1          ggridges_0.5.2         fastmap_1.1.0          yaml_2.2.1             goftest_1.2-2          bit64_4.0.5            fitdistrplus_1.1-1    
# [81] purrr_0.3.4            RANN_2.6.1             pbapply_1.4-3          future_1.19.1          nlme_3.1-149           mime_0.9               grr_0.9.5              compiler_4.0.2        
# [89] rstudioapi_0.11        plotly_4.9.2.1         png_0.1-7              spatstat.utils_2.2-0   tibble_3.0.4           statmod_1.4.35         geneplotter_1.66.0     stringi_1.7.6         
# [97] forcats_0.5.0          lattice_0.20-41        nloptr_1.2.2.2         vctrs_0.4.1            pillar_1.4.6           lifecycle_0.2.0        spatstat.geom_2.4-0    lmtest_0.9-38         
# [105] RcppAnnoy_0.0.19       data.table_1.13.0      cowplot_1.1.0          bitops_1.0-6           irlba_2.3.3            httpuv_1.5.4           patchwork_1.0.1        R6_2.4.1              
# [113] promises_1.1.1         KernSmooth_2.23-17     gridExtra_2.3          codetools_0.2-16       boot_1.3-25            MASS_7.3-53            withr_2.3.0            sctransform_0.3.2     
# [121] GenomeInfoDbData_1.2.3 mgcv_1.8-33            grid_4.0.2             rpart_4.1-15           glmmTMB_1.1.3          tidyr_1.1.2            minqa_1.2.4            Rtsne_0.15            
# [129] numDeriv_2016.8-1.1    shiny_1.5.0   
