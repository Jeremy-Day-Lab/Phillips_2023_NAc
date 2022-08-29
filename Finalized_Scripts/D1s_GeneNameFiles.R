#Create files for HOMER findmotifs.pl
setwd("/data/project/daylab/2019-JD-0040/MCN_Code/Tables/DESeq2_D1s/")
set.seed(1234)

#Drd1-MSN-2
D1_2_enriched <- read.table("Drd1_DESeq2_Down.txt",header = TRUE)
write.table(x         = data.frame(GeneName = D1_2_enriched$GeneName),
            file      = "Drd1_MSN_2_enriched_genesonly.txt",
            sep       = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE)

#Drd1-MSN-1
D1_1_enriched <- read.table("Drd1_DESeq2_Up.txt",header = TRUE)
write.table(x         = data.frame(GeneName = D1_1_enriched$GeneName),
            file      = "Drd1_MSN_1_enriched_genesonly.txt",
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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# loaded via a namespace (and not attached):
#   [1] compiler_4.0.2 tools_4.0.2    yaml_2.2.1    
