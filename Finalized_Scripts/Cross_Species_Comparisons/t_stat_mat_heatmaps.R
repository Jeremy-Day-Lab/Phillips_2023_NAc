#Set working directory and seed
setwd("/data/project/daylab/2019-JD-0040/MCN_Code")
set.seed(1234)

######Load libraries
library(Seurat) #For the NoGrid() function. 
library(ggplot2)

#Read in the rat t stat mat 
Rat_t_stat_mat <- read.table(file = "Rn7_Neuron_DEG_Orthologs_t_stat_mat.txt",sep = "\t")
Rat_t_stat_mat <- Rat_t_stat_mat[c("Drd1","Ebf1","Htr4","Drd2","Penk","Drd3","Grm8","Oprm1"),
                                 c("Drd1.MSN.1","Drd1.MSN.2","Drd2.MSN.1","Drd2.MSN.2","Drd3.MSN","Grm8.MSN")]
Rat_t_stat_mat$Gene <- row.names(Rat_t_stat_mat)
Rat_t_melt <- reshape2::melt(Rat_t_stat_mat)
Rat_t_melt$value <- ifelse(Rat_t_melt$value >= 50,
                           50,
                           Rat_t_melt$value)
Rat_t_melt$value <- ifelse(Rat_t_melt$value <= (-50),
                           -50,
                           Rat_t_melt$value)
Rat_t_melt$Gene <- factor(x = Rat_t_melt$Gene,
                          levels = rev(c("Drd1","Ebf1","Htr4","Drd2","Penk","Drd3","Grm8","Oprm1")))
Rat <- ggplot(data = Rat_t_melt,aes(x = variable,y= Gene,fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = c("dodgerblue","lightgray","red")) +
  theme_bw() +
  ggtitle("Rat MSNs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(hjust = 0.5)) +
  NoGrid() +
  labs(x    = "Cell Type",
       fill = "t-statistic")

ggsave(filename = "Plots/Rat_t_stat_marker_gene_heatmap.pdf",
       plot     = Rat,
       height   = 8,
       width    = 8)

#Read in the NHP t stat mat 
NHP_t_stat_mat <- read.table(file = "NHP_MSNs_DEG_Orthologs_t_stat_mat.txt",sep = "\t")
NHP_t_stat_mat <- NHP_t_stat_mat[c("DRD1","EBF1","HTR4","DRD2","PENK","DRD3","GRM8","OPRM1"),]
NHP_t_stat_mat$Gene <- row.names(NHP_t_stat_mat)
NHP_t_melt <- reshape2::melt(NHP_t_stat_mat)
NHP_t_melt$value <- ifelse(NHP_t_melt$value >= 20,
                           20,
                           NHP_t_melt$value)
NHP_t_melt$value <- ifelse(NHP_t_melt$value <= (-20),
                           -20,
                           NHP_t_melt$value)
NHP_t_melt$Gene <- factor(x = NHP_t_melt$Gene,
                          levels = rev(c("DRD1","EBF1","HTR4","DRD2","PENK","DRD3","GRM8","OPRM1")))
NHP <- ggplot(data = NHP_t_melt,aes(x = variable,y= Gene,fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = c("dodgerblue","lightgray","red")) +
  theme_bw() +
  ggtitle("NHP MSNs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(hjust = 0.5)) +
  NoGrid() +
  labs(x    = "Cell Type",
       fill = "t-statistic")

ggsave(filename = "Plots/NHP_t_stat_marker_gene_heatmap.pdf",
       plot     = NHP,
       height   = 8,
       width    = 8)

#Read in the human t stat mat 
human_t_stat_mat <- read.table(file = "human_neurons_DEG_Orthologs_t_stat_mat.txt",sep = "\t")
human_t_stat_mat <- human_t_stat_mat[c("DRD1","EBF1","HTR4","DRD2","PENK","DRD3","GRM8","OPRM1"),
                                     c("MSN.D1_A","MSN.D1_B","MSN.D1_C","MSN.D1_D","MSN.D1_E","MSN.D1_F", 
                                       "MSN.D2_A","MSN.D2_B","MSN.D2_C","MSN.D2_D")]
human_t_stat_mat$Gene <- row.names(human_t_stat_mat)
human_t_melt <- reshape2::melt(human_t_stat_mat)
human_t_melt$value <- ifelse(human_t_melt$value >= 50,
                             50,
                             human_t_melt$value)
human_t_melt$value <- ifelse(human_t_melt$value <= (-50),
                             -50,
                             human_t_melt$value)
human_t_melt$Gene <- factor(x = human_t_melt$Gene,
                            levels = rev(c("DRD1","EBF1","HTR4","DRD2","PENK","DRD3","GRM8","OPRM1")))
human <- ggplot(data = human_t_melt,aes(x = variable,y= Gene,fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = c("dodgerblue","lightgray","red")) +
  theme_bw() +
  ggtitle("Human MSNs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(hjust = 0.5)) +
  NoGrid() +
  labs(x    = "Cell Type",
       fill = "t-statistic")

ggsave(filename = "Plots/Human_t_stat_marker_gene_heatmap.pdf",
       plot     = human,
       height   = 8,
       width    = 8)

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
#   [1] nlme_3.1-149          matrixStats_0.57.0    spatstat.sparse_2.0-0 RcppAnnoy_0.0.19      RColorBrewer_1.1-2   
# [6] httr_1.4.2            sctransform_0.3.2     tools_4.0.2           R6_2.4.1              irlba_2.3.3          
# [11] rpart_4.1-15          KernSmooth_2.23-17    uwot_0.1.10           mgcv_1.8-33           lazyeval_0.2.2       
# [16] colorspace_1.4-1      withr_2.3.0           tidyselect_1.1.0      gridExtra_2.3         compiler_4.0.2       
# [21] cli_3.3.0             plotly_4.9.2.1        labeling_0.3          scales_1.1.1          lmtest_0.9-38        
# [26] spatstat.data_2.1-0   ggridges_0.5.2        pbapply_1.4-3         rappdirs_0.3.1        goftest_1.2-2        
# [31] stringr_1.4.0         digest_0.6.26         spatstat.utils_2.2-0  pkgconfig_2.0.3       htmltools_0.5.2      
# [36] fastmap_1.1.0         htmlwidgets_1.5.2     rlang_1.0.2           rstudioapi_0.11       shiny_1.5.0          
# [41] farver_2.0.3          generics_0.0.2        zoo_1.8-8             jsonlite_1.7.1        ica_1.0-2            
# [46] dplyr_1.0.2           magrittr_1.5          patchwork_1.0.1       Matrix_1.3-4          Rcpp_1.0.7           
# [51] munsell_0.5.0         abind_1.4-5           reticulate_1.16       lifecycle_0.2.0       stringi_1.7.6        
# [56] yaml_2.2.1            MASS_7.3-53           Rtsne_0.15            plyr_1.8.6            grid_4.0.2           
# [61] parallel_4.0.2        listenv_0.8.0         promises_1.1.1        ggrepel_0.8.2         crayon_1.3.4         
# [66] miniUI_0.1.1.1        deldir_1.0-6          lattice_0.20-41       cowplot_1.1.0         splines_4.0.2        
# [71] tensor_1.5            pillar_1.4.6          igraph_1.2.6          spatstat.geom_2.4-0   future.apply_1.6.0   
# [76] reshape2_1.4.4        codetools_0.2-16      leiden_0.3.3          glue_1.6.2            data.table_1.13.0    
# [81] png_0.1-7             vctrs_0.4.1           httpuv_1.5.4          gtable_0.3.0          RANN_2.6.1           
# [86] purrr_0.3.4           spatstat.core_2.3-0   polyclip_1.10-0       tidyr_1.1.2           scattermore_0.7      
# [91] future_1.19.1         mime_0.9              xtable_1.8-4          later_1.1.0.1         survival_3.2-7       
# [96] viridisLite_0.3.0     tibble_3.0.4          cluster_2.1.0         globals_0.13.1        fitdistrplus_1.1-1   
# [101] ellipsis_0.3.2        ROCR_1.0-11          
