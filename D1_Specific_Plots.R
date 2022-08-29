#Set working directory and seed
setwd("/data/project/daylab/")
set.seed(1234)

######Load libraries
library(Seurat)
library(ggplot2)
library(scales)

#Load the NAc_Combo object
NAc_Combo <- readRDS(file = "2019-JD-0040/MCN_Code/Objects/NAc_Combo_Integrated.RDS")

DimPlot(object = NAc_Combo,reduction = "umap",label = TRUE) + NoLegend()


#Read in the D1_1 and D1_2 enriched files
D1_1_enriched <- read.table("2019-JD-0040/MCN_Code/Tables/DESeq2_D1s/Drd1_DESeq2_Up.txt",header = TRUE)
D1_2_enriched <-  read.table("2019-JD-0040/MCN_Code/Tables/DESeq2_D1s/Drd1_DESeq2_Down.txt",header = TRUE)

#All Pct_Change in D1_2_enriched pops is negative due to calculation previously used. Need to make absolute value. 
D1_2_enriched$Abs_Pct_Change <- abs(D1_2_enriched$Pct_Change)

#Get Average expression values for every gene
D1s <- subset(NAc_Combo,idents = c("Drd1-MSN-1","Drd1-MSN-2"))
DefaultAssay(D1s) <- "RNA"
D1_Averages <- as.data.frame(AverageExpression(object = D1s,assays = "RNA",slot = "data")$RNA)
D1_Averages$Gene <- row.names(D1s)

#Merge
D1_2_enriched <- merge(x = D1_2_enriched,
                       y = D1_Averages[,c("Drd1-MSN-2","Gene")],
                       by.x = "GeneName",
                       by.y = "Gene")
colnames(D1_2_enriched)[12] <- "CountValues"


D1_Enriched <- read.table("2019-JD-0040/MCN_Code/Tables/DESeq2_D1s/Drd1_DESeq2_results.txt",header = TRUE)
D1_Enriched <- D1_Enriched[order(D1_Enriched$Pct_Change,decreasing = FALSE),]
D1_Enriched$Rank <- 1:nrow(D1_Enriched)
D1_Enriched <-  merge(x = D1_Enriched,
                      y = D1_Averages,
                      by.x = "GeneName",
                      by.y = "Gene")
D1_Enriched$Max_Expression <- apply(X = D1_Enriched[,c("Drd1-MSN-1","Drd1-MSN-2")],MARGIN = 1,FUN = max)

#Create a named that includes the number 
colors <- hue_pal()(16)
names(colors) <- levels(NAc_Combo$Combo_CellType)

#Read in  the dataframe containing gini calculations
Gini_df <- read.table(file = "2019-JD-0040/MCN_Code/Tables/Gini_Dataframe_AllClusters.txt",sep = "\t",header = TRUE)


D1_2_Gini <- merge(x = D1_Enriched,
                   y = Gini_df,
                   by.x = "GeneName",
                   by.y = "Gene")
D1_2_Gini$logCluster <- log10(D1_2_Gini[,"Drd1.MSN.2"])
D1_2_Gini <- D1_2_Gini[!is.infinite(D1_2_Gini$logCluster),]

#Get percent calculation
cells.1    <- WhichCells(NAc_Combo,idents = "Drd1-MSN-2")
counts_mat <- NAc_Combo@assays$RNA@counts
D1_2_Gini$Pct_Expressing <- NA
D1_2_Gini$Pct_Expressing <- (rowSums(x = counts_mat[D1_2_Gini$GeneName, cells.1,drop = FALSE] > 0) / length(cells.1))*100

#MAke a gini plot using the Gini calcs from all clusters and D1/D2 DEGs
D1_2_Gini_plot <- ggplot(data = D1_2_Gini,aes(x = logCluster,y = Gini_calc,size = Pct_Expressing)) +
  geom_point(color = "lightgrey",aes(size = Pct_Expressing)) +
  geom_point(data = subset(D1_2_Gini,subset=(padj <= 0.05 & log2FoldChange < 0)),
             color = "dodgerblue") +
  geom_point(data = D1_2_Gini[D1_2_Gini$GeneName %in% c("Slit2","Cyct","Nell1","Trhde","Col11a1","Htr4","Oprd1","Trpc7","Kcnt2","AABR07007642.1"),],
             color = "limegreen") +
  geom_text_repel(data  = D1_2_Gini[D1_2_Gini$GeneName %in% c("Slit2","Cyct","Nell1","Trhde","Col11a1","Htr4","Oprd1","Trpc7","Kcnt2","AABR07007642.1"),],
                  label = D1_2_Gini[D1_2_Gini$GeneName %in% c("Slit2","Cyct","Nell1","Trhde","Col11a1","Htr4","Oprd1","Trpc7","Kcnt2","AABR07007642.1"),"GeneName"],
                  show.legend = FALSE,
                  size        = 6) +
  theme_bw() +
  NoGrid() + 
  labs(x     = "log10(Drd1-MSN-2 Counts)",
       y     = "Gini Coefficient",
       size  = "% of Cells in Cluster\nExpressing Gene",
       title = "Drd1-MSN-2") +
  theme(plot.title = element_text(hjust = 0.5))

#save the plot. 
ggsave(plot     = D1_2_Gini_plot,
       filename = "2019-JD-0040/MCN_Code/Plots/Drd1_MSN_2_GiniCalcs.pdf",
       height   = 12,
       width    = 12)
  



#Order the dataframe by the rank of the 
library(ggrepel)
ggplot() +
  geom_point(data = subset(D1_Enriched,subset=(log2FoldChange > 0)),aes(x = log2FoldChange, y = Pct_Change,size = D1_1_Pct_Expressing),color = colors["Drd1-MSN-1"]) +
  geom_point(data = subset(D1_Enriched,subset=(log2FoldChange < 0)),aes(x = log2FoldChange, y = Pct_Change,size = D1_2_Pct_Expressing),color = colors["Drd1-MSN-2"]) +
  xlim(c(-6,6)) +
  ylim(c(-100,100)) +
  theme_bw() +
  NoGrid() +
  labs(x = "log2(Drd1-MSN-1/Drd1-MSN-2)",
       y = "% of Drd1-MSN-1 expressing - % of Drd1-MSN-2 expressing",
       size = "% of Cells \nin cluster expressing") +
  geom_hline(yintercept = 0,lty = 2) +
  geom_vline(xintercept = 0,lty = 2) +
  geom_text_repel(data = subset(D1_Enriched,subset=(log2FoldChange > 0))[order(subset(D1_Enriched,subset=(log2FoldChange > 0))$Pct_Change,decreasing = TRUE)[1:10],],
                  aes(x = log2FoldChange, y = Pct_Change),
                  label = subset(D1_Enriched,subset=(log2FoldChange > 0))[order(subset(D1_Enriched,subset=(log2FoldChange > 0))$Pct_Change,decreasing = TRUE)[1:10],"GeneName"]) +
  geom_text_repel(data = subset(D1_Enriched,subset=(log2FoldChange < 0))[order(subset(D1_Enriched,subset=(log2FoldChange < 0))$Pct_Change,decreasing = FALSE)[1:10],],
                  aes(x = log2FoldChange, y = Pct_Change),
                  label = subset(D1_Enriched,subset=(log2FoldChange < 0))[order(subset(D1_Enriched,subset=(log2FoldChange < 0))$Pct_Change,decreasing = FALSE)[1:10],"GeneName"])

FeaturePlot(object = NAc_Combo,features = c("Drd1","Ebf1","Calb1","Htr4","Chst9"),cols = c("lightgrey","red"),ncol = 2,order = TRUE)
VlnPlot(object = NAc_Combo,features = c("Drd1","Ebf1","Calb1","Htr4","Reln"),stack = TRUE,flip = TRUE)



#New subset
Cluster_Lists <- vector(mode = "list",length = length(levels(NAc_Combo)))
#name every element of the list the numeric identity from subclustering
names(Cluster_Lists) <- levels(NAc_Combo)
#Within each element of the list, create a dataframe that will contain the percentage of cells within the cluster expressing two genes 
for(i in names(Cluster_Lists)){
  Cluster_Lists[[i]] <- data.frame(Drd1      = rep(NA,13),#
                                   Ebf1      = rep(NA,13),#
                                   Calb1     = rep(NA,13),#
                                   Pdyn      = rep(NA,13),#
                                   Reln      = rep(NA,13),#
                                   Htr4      = rep(NA,13),#
                                   Nell1     = rep(NA,13),#
                                   Slit2     = rep(NA,13),#
                                   Drd2      = rep(NA,13),#
                                   Penk      = rep(NA,13),#
                                   Drd3      = rep(NA,13),#
                                   Grm8      = rep(NA,13),#
                                   Chst9     = rep(NA,13),#
                                   row.names = c("Drd1","Ebf1","Calb1","Pdyn","Reln","Htr4","Nell1",
                                                 "Slit2","Drd2","Penk","Drd3","Grm8","Chst9"))
  
}

#Calculate the percentage of cells expressing two genes 
for(i in levels(NAc_Combo)){
  print(i)
  #Subset the VTA_subset object for a single identity
  NAc_Combo2 <- subset(NAc_Combo,idents = i)
  #Pull the count values for the 13 genes and transpose the dataframe so that cells are rows and genes are columns 
  x <- as.data.frame(t(as.matrix(GetAssayData(object = NAc_Combo2,assay = "RNA")[c("Drd1","Ebf1","Calb1","Pdyn","Reln","Htr4","Nell1",
                                                                                   "Slit2","Drd2","Penk","Drd3","Grm8","Chst9"),])))
  #Now loop through the 13 genes of interest by looping through the column names of each dataframe created in the first step 
  for(l in colnames(Cluster_Lists[[as.character(i)]])){
    print(l)
    #Get all of the genes that are not equal to the iteration of the loop
    Other_Genes <- colnames(Cluster_Lists[[as.character(i)]])[which(colnames(Cluster_Lists[[as.character(i)]])!=l)]
    for(k in Other_Genes){
      #Calculate the percentage of cells expressing the l and k, both of which are some combination of the 13 genes
      #Dividing by ncol(NAc_Combo2) divides by the number of cells within the cluster of interest 
      Cluster_Lists[[as.character(i)]][l,k] <- (length(which(x[,l] > 0 & x[,k] > 0))/ncol(NAc_Combo2))*100
    }
  }
}


library(pheatmap)
for(i in names(Cluster_Lists)){
  pdf(file = paste0("/data/project/daylab/2019-JD-0040/MCN_Code/Plots/CoexpressionHeatmaps/",i,"_heatmap.pdf"),
      height = 8,
      width = 8)
  breaksList <- seq(0, 100, by = 5)
  pheatmap(mat = Cluster_Lists[[i]],
           color = colorRampPalette(c("gray96","peachpuff","goldenrod1","orange","darkorange","red","red4"))(length(breaksList)),
           breaks = breaksList,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           main = i)
  dev.off()
}


  
  
