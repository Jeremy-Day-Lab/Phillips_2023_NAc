#There are some interesting differences between the Drd1 populations. Mainly that Drd1-MSN-2 is devoid of Ebf1. 
#Goal of this analysis is to find a gene specific to Drd1.MSN.2
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

#Subset NAc for only Drd1
D1s <- subset(NAc_Combo, subset=(Combo_CellType == "Drd1-MSN-1" |
                                   Combo_CellType == "Drd1-MSN-2"))

DimPlot(object = D1s,reduction = "umap",label = TRUE) + NoLegend()


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
D1s$sample.id <- as.factor(paste(D1s$Dataset,D1s$Sex_Stim,sep = "_"))

#Get cell and sample metrics for aggregation
D1s@meta.data$Combo_CellType <- as.character(D1s@meta.data$Combo_CellType)
D1s@meta.data$CellType <- as.factor(D1s@meta.data$Combo_CellType)
cell_names <- purrr::set_names(levels(D1s@meta.data$CellType))
cell_names
# Drd1-MSN-1   Drd1-MSN-2 
# "Drd1-MSN-1" "Drd1-MSN-2" 

#Number of clusters
cluster_total <- length(cell_names)
cluster_total
# [1] 2

# Named vector of sample names
D1s$Sex_Stim <- factor(D1s$Sex_Stim)
sample_ids <- purrr::set_names(levels(D1s@meta.data$Sex_Stim))
sample_ids


# Total number of samples 
sample_total <- length(sample_ids)
sample_total
# Fem_Coc    Fem_Sal   Male_Coc   Male_Sal 
# "Fem_Coc"  "Fem_Sal" "Male_Coc" "Male_Sal" 

#Figure out how many cells in each main group
table(D1s@meta.data$Sex_Stim)
# Fem_Coc  Fem_Sal Male_Coc Male_Sal 
# 2259     1482     1887     1674 

##########Count aggregation to sample level
groups <- D1s@meta.data[, c("CellType", "sample.id")]

# Aggregate across cluster-sample groups (raw counts, thus counts slot)
count_aggr <- Matrix.utils::aggregate.Matrix(t(D1s@assays$RNA@counts), 
                                             groupings = groups, fun = "sum") 

dim(count_aggr)
# [1]    16 30560

#Transpose count_aggr
count_aggr_t <- t(count_aggr)

#change the count_aggr_t
colnames(count_aggr_t) <- gsub(x = colnames(count_aggr_t),pattern = "-",replacement = ".")


# Create a data frame with the sample IDs, cluster IDs and condition
metadata <- data.frame(cluster_id = sub("_.*","", colnames(count_aggr_t)),
                       sample_id = colnames(count_aggr_t),
                       Sex.Stim = sub("^[^_]*_","", colnames(count_aggr_t)))
#Create extra columns in the metadata
metadata$dataset <- as.factor(as.character(lapply(strsplit(metadata$Sex.Stim,"_"),"[",1)))
metadata$Sex     <- as.factor(as.character(lapply(strsplit(metadata$Sex.Stim,"_"),"[",2)))
metadata$Stim    <- as.factor(as.character(lapply(strsplit(metadata$Sex.Stim,"_"),"[",3)))


#Create the DESeq2 object
dds <- DESeqDataSetFromMatrix(count_aggr_t, 
                              colData = metadata, 
                              design = ~ dataset + Stim + cluster_id)
dds$cluster_id <- relevel(dds$cluster_id, ref = "Drd1.MSN.2")
keep <- rowMeans(counts(dds)) > 5
dds <- dds[keep,]

############Quality control############
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

#Create a dataframe that includes PC information
data <- plotPCA(vsd, intgroup = c("cluster_id"), returnData = TRUE)
#Calculate percent variance
percentVar <- round(100 * attr(data, "percentVar"))
#Create columns for dataset, sex, and stim
data$dataset <- as.character(lapply(strsplit(x = data$name,split = "_"),"[",2))
data$Sex     <- as.character(lapply(strsplit(x = data$name,split = "_"),"[",3))
data$Stim    <- as.character(lapply(strsplit(x = data$name,split = "_"),"[",4))
#Create a sex_stim column
data$Sex_Stim <- paste(data$Sex,data$Stim,sep = "_")

#PCA for major celltype
pca.plot.CellType <- ggplot(data, aes(PC1, PC2, color = cluster_id)) +
  geom_point(size = 7) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme(text = element_text(size = 20)) +
  theme_bw(base_size = 16)

#PCA for dataet
pca.plot.dataset <- ggplot(data, aes(PC1, PC2, color = dataset)) +
  geom_point(size = 7) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme(text = element_text(size = 20)) +
  theme_bw(base_size = 16)

#PCA for stim
pca.plot.stim <- ggplot(data, aes(PC1, PC2, color = Stim)) +
  geom_point(size = 7) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme(text = element_text(size = 20)) +
  theme_bw(base_size = 16)

#PCA for sex
pca.plot.sex <- ggplot(data, aes(PC1, PC2, color = Sex)) +
  geom_point(size = 7) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme(text = element_text(size = 20)) +
  theme_bw(base_size = 16)

#PCA for sex.stim
pca.plot.sex_stim <- ggplot(data, aes(PC1, PC2, color = Sex_Stim)) +
  geom_point(size = 7) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme(text = element_text(size = 20)) +
  theme_bw(base_size = 16)

ggsave(plot = cowplot::plot_grid(plotlist = list(pca.plot.CellType,
                                                 pca.plot.dataset,
                                                 pca.plot.stim,
                                                 pca.plot.sex,
                                                 pca.plot.sex_stim),
                                 ncol = 3),
       filename = "2019-JD-0040/MCN_Code/Plots/DESeq2_D1s/PCA_Plots_Combo_withStimasCovariate.pdf",
       height = 8,
       width = 16)

#Now run DESeq2 workflow
dds <- DESeq(dds, test="LRT", reduced = ~ dataset + Stim)
res <- results(dds)
res_df <- as.data.frame(res)
res_df$GeneName <- rownames(res_df)


#Calculate the percent of cells in each cluster that express these genes
#Drd1-MSN-1
cells.1    <- WhichCells(D1s,idents = "Drd1-MSN-1")
counts_mat <- D1s@assays$RNA@counts
res_df$D1_1_Pct_Expressing <- NA
res_df$D1_1_Pct_Expressing <- (rowSums(x = counts_mat[res_df$GeneName, cells.1,drop = FALSE] > 0) / length(cells.1))*100

#Drd1-MSN-2
cells.2    <- WhichCells(D1s,idents = "Drd1-MSN-2")
res_df$D1_2_Pct_Expressing <- NA
res_df$D1_2_Pct_Expressing <- (rowSums(x = counts_mat[res_df$GeneName, cells.2,drop = FALSE] > 0) / length(cells.2))*100


#Calculate difference in percent expressing
res_df$Pct_Change <- res_df$D1_1_Pct_Expressing - res_df$D1_2_Pct_Expressing

#Make a volcano plot
library(ggrepel)
VP1 <- ggplot(data = res_df,aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point() +
  geom_point(data = subset(res_df,subset=(padj <= 0.05 & log2FoldChange >= 0.5)),color = "tomato") +
  geom_point(data = subset(res_df,subset=(padj <= 0.05 & log2FoldChange <= (-0.5))),color = "dodgerblue") +
  xlim(c(-6,6)) +
  ylim(c(0,220)) +
  geom_hline(yintercept = -log10(0.05),lty = 2) +
  geom_vline(xintercept = 0.5,lty = 2) +
  geom_vline(xintercept = -0.5,lty = 2) +
  geom_text_repel(data  = subset(res_df,subset=(padj <= 0.05 & log2FoldChange <0))[order(subset(res_df,subset=(padj <= 0.05 & log2FoldChange <0))$padj,decreasing=FALSE)[1:10],],
                  label = subset(res_df,subset=(padj <= 0.05 & log2FoldChange <0))[order(subset(res_df,subset=(padj <= 0.05 & log2FoldChange <0))$padj,decreasing=FALSE)[1:10],"GeneName"]) +
  geom_text_repel(data  = subset(res_df,subset=(padj <= 0.05 & log2FoldChange >0))[order(subset(res_df,subset=(padj <= 0.05 & log2FoldChange >0))$padj,decreasing=FALSE)[1:10],],
                  label = subset(res_df,subset=(padj <= 0.05 & log2FoldChange >0))[order(subset(res_df,subset=(padj <= 0.05 & log2FoldChange >0))$padj,decreasing=FALSE)[1:10],"GeneName"]) +
  geom_text_repel(data  = subset(res_df,subset=(GeneName == "Ebf1")),
                  label = subset(res_df,subset=(GeneName == "Ebf1"))$GeneName) +
  theme_bw() +
  NoGrid() +
  annotate("text",x = -3.5,y = 215,label = paste("Enriched in Drd1-MSN-2\n",nrow(subset(res_df,subset=(padj <= 0.05 & log2FoldChange <(-0.5)))),"genes"),color = "dodgerblue") +
  annotate("text",x = 3.5,y = 215,label = paste("Enriched in Drd1-MSN-1\n",nrow(subset(res_df,subset=(padj <= 0.05 & log2FoldChange >0.5))),"genes"),color = "tomato") +
  labs(x = "log2(Drd1-MSN-1/Drd1-MSN-2)",
       y = "-log10(Adjusted p-value)")
ggsave(plot     = VP1,
       filename = "2019-JD-0040/MCN_Code/Plots/DESeq2_D1s/VolcanoPlot_StimasCovariate.pdf",
       height   = 8,
       width    = 8)

#Calculate Percent expressing
res_df$Pct_Change <- res_df$D1_1_Pct_Expressing - res_df$D1_2_Pct_Expressing
Fp1 <- FeaturePlot(object = NAc_Combo,features = c("Drd1","Ebf1","Calb1","Htr4","Camk2d","Nell1"),cols = c("lightgrey","red"),ncol =2)
ggsave(plot     = Fp1,
       filename = "2019-JD-0040/MCN_Code/Plots/DESeq2_D1s/TopGenes_FeaturePlot.pdf",
       height   = 12,
       width    = 8)


#Write out the tables before going any further. 
write.table(x         = res_df,
            file      = "2019-JD-0040/MCN_Code/Tables/DESeq2_D1s/Drd1_DESeq2_results_Stimascovariate.txt",
            sep       = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE)

write.table(x         = subset(res_df,subset=(padj <= 0.05 & log2FoldChange < (-0.5))),
            file      = "2019-JD-0040/MCN_Code/Tables/DESeq2_D1s/Drd1_DESeq2_Down_Stimascovariate.txt",
            sep       = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE)

write.table(x         = subset(res_df,subset=(padj <= 0.05 & log2FoldChange > 0.5)),
            file      = "2019-JD-0040/MCN_Code/Tables/DESeq2_D1s/Drd1_DESeq2_Up_Stimascovariate.txt",
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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggrepel_0.8.2               pheatmap_1.0.12             RColorBrewer_1.1-2          DESeq2_1.28.1               SummarizedExperiment_1.18.2
# [6] DelayedArray_0.14.1         matrixStats_0.57.0          Biobase_2.48.0              GenomicRanges_1.40.0        GenomeInfoDb_1.24.2        
# [11] IRanges_2.22.2              S4Vectors_0.26.1            BiocGenerics_0.34.0         Matrix.utils_0.9.8          Matrix_1.3-4               
# [16] ComplexUpset_1.3.3          stringr_1.4.0               dplyr_1.0.2                 Libra_1.0.0                 ggplot2_3.3.2              
# [21] SeuratObject_4.0.2          Seurat_4.0.4               
# 
# loaded via a namespace (and not attached):
#   [1] blme_1.0-5             plyr_1.8.6             igraph_1.2.6           lazyeval_0.2.2         TMB_1.8.1              splines_4.0.2         
# [7] BiocParallel_1.22.0    listenv_0.8.0          scattermore_0.7        digest_0.6.26          htmltools_0.5.2        lmerTest_3.1-3        
# [13] magrittr_1.5           memoise_1.1.0          tensor_1.5             cluster_2.1.0          ROCR_1.0-11            limma_3.44.3          
# [19] globals_0.13.1         annotate_1.66.0        tester_0.1.7           spatstat.sparse_2.0-0  colorspace_1.4-1       blob_1.2.1            
# [25] rappdirs_0.3.1         crayon_1.3.4           RCurl_1.98-1.2         jsonlite_1.7.1         genefilter_1.70.0      lme4_1.1-26           
# [31] spatstat.data_2.1-0    survival_3.2-7         zoo_1.8-8              glue_1.6.2             polyclip_1.10-0        gtable_0.3.0          
# [37] zlibbioc_1.34.0        XVector_0.28.0         leiden_0.3.3           future.apply_1.6.0     abind_1.4-5            scales_1.1.1          
# [43] DBI_1.1.0              edgeR_3.30.3           miniUI_0.1.1.1         Rcpp_1.0.7             viridisLite_0.3.0      xtable_1.8-4          
# [49] reticulate_1.16        spatstat.core_2.3-0    bit_4.0.4              htmlwidgets_1.5.2      httr_1.4.2             ellipsis_0.3.2        
# [55] ica_1.0-2              farver_2.0.3           pkgconfig_2.0.3        XML_3.99-0.5           uwot_0.1.10            deldir_1.0-6          
# [61] locfit_1.5-9.4         labeling_0.3           tidyselect_1.1.0       rlang_1.0.2            reshape2_1.4.4         later_1.1.0.1         
# [67] AnnotationDbi_1.50.3   pbmcapply_1.5.0        munsell_0.5.0          tools_4.0.2            cli_3.3.0              generics_0.0.2        
# [73] RSQLite_2.2.1          ggridges_0.5.2         fastmap_1.1.0          yaml_2.2.1             goftest_1.2-2          bit64_4.0.5           
# [79] fitdistrplus_1.1-1     purrr_0.3.4            RANN_2.6.1             pbapply_1.4-3          future_1.19.1          nlme_3.1-149          
# [85] mime_0.9               grr_0.9.5              compiler_4.0.2         rstudioapi_0.11        plotly_4.9.2.1         png_0.1-7             
# [91] spatstat.utils_2.2-0   tibble_3.0.4           statmod_1.4.35         geneplotter_1.66.0     stringi_1.7.6          forcats_0.5.0         
# [97] lattice_0.20-41        nloptr_1.2.2.2         vctrs_0.4.1            pillar_1.4.6           lifecycle_0.2.0        spatstat.geom_2.4-0   
# [103] lmtest_0.9-38          RcppAnnoy_0.0.19       data.table_1.13.0      cowplot_1.1.0          bitops_1.0-6           irlba_2.3.3           
# [109] httpuv_1.5.4           patchwork_1.0.1        R6_2.4.1               promises_1.1.1         KernSmooth_2.23-17     gridExtra_2.3         
# [115] codetools_0.2-16       boot_1.3-25            MASS_7.3-53            withr_2.3.0            sctransform_0.3.2      GenomeInfoDbData_1.2.3
# [121] mgcv_1.8-33            grid_4.0.2             rpart_4.1-15           glmmTMB_1.1.3          tidyr_1.1.2            minqa_1.2.4           
# [127] Rtsne_0.15             numDeriv_2016.8-1.1    shiny_1.5.0           
