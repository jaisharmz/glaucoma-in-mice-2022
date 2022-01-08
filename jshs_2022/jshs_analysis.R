# Import necessary libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(readxl)
library(cellranger)
library(SingleCellExperiment)
library(TSCAN)
library(scater)
library(ggplot2)
library(ggthemes)
library(ggbeeswarm)
library(corrplot)
library(scuttle)
library(gplots)
library(RColorBrewer)
library(rhdf5)
library(mclust)
library(igraph)
library(wesanderson)
library(destiny)
library(slingshot)
library(monocle)
library(edgeR)
library(DESeq2)
library("xlsx")

# Set working directory
setwd("/Users/jaisharma/Documents/jshs_2022")

# Define important functions
analysis <- function(group1, group2){
  super_directory <- paste(group1,"_vs_",group2,sep="")
  dir.create(super_directory)
  
  cf_et <- exactTest(cf_all, pair = c(group1, group2))
  DE <- topTags(cf_et , n = nrow(cf_et$table))$table
  
  # DE[match(c("C1qa", "C1qb", "C1qc"), DE$hgnc_symbol),]
  
  # Upregulated genes
  up_regulated <- DE[DE$logFC>0,]
  up_regulated <- up_regulated[order(up_regulated$logFC, decreasing=T),c(2,3,5,6)]
  up_regulated_20 <- up_regulated[c(1:20),]
  up_regulated_50 <- up_regulated[c(1:50),]
  write.csv(up_regulated, paste(super_directory, "/", group1,"_vs_",group2,"_de_all_up.csv",sep=""), row.names=T)
  write.csv(up_regulated_20, paste(super_directory, "/", group1,"_vs_",group2,"_de20_up.csv",sep=""), row.names=T)
  write.csv(up_regulated_50, paste(super_directory, "/", group1,"_vs_",group2,"_de50_up.csv",sep=""), row.names=T)
  
  # Downregulated genes
  down_regulated <- DE[DE$logFC<0,]
  down_regulated <- down_regulated[order(down_regulated$logFC, decreasing=F),c(2,3,5,6)]
  down_regulated_20 <- down_regulated[c(1:20),]
  down_regulated_50 <- down_regulated[c(1:50),]
  write.csv(down_regulated, paste(super_directory, "/", group1,"_vs_",group2,"_de_all_down.csv",sep=""), row.names=T)
  write.csv(down_regulated_20, paste(super_directory, "/", group1,"_vs_",group2,"_de20_down.csv",sep=""), row.names=T)
  write.csv(down_regulated_50, paste(super_directory, "/", group1,"_vs_",group2,"_de50_down.csv",sep=""), row.names=T)
  
  # Heatmap setup
  genes_subset <- c(up_regulated_20$hgnc_symbol[c(1:10)], down_regulated_20$hgnc_symbol[c(1:10)])
  genes_subset_index <- match(genes_subset, genes_all)
  col_index_group1 <- group == group1
  matrix_nums_heatmap_group1 <- cpmTMM[match(genes_subset, genes_all), c(1:(dim(cpmTMM)[2]-2))][,col_index_group1]
  matrix_nums_heatmap_group1 <- matrix_nums_heatmap_group1[,order(colSums(matrix_nums_heatmap_group1))]
  col_index_group2 <- group == group2
  matrix_nums_heatmap_group2 <- cpmTMM[match(genes_subset, genes_all), c(1:(dim(cpmTMM)[2]-2))][,col_index_group2]
  matrix_nums_heatmap_group2 <- matrix_nums_heatmap_group2[,order(colSums(matrix_nums_heatmap_group2))]
  matrix_nums_heatmap <- as.matrix(cbind(matrix_nums_heatmap_group1,matrix_nums_heatmap_group2))
  matrix_gene_heatmap <- as.matrix(cpmTMM$hgnc_symbol[genes_subset_index])
  title <- paste(group1," vs ",group2,sep="")
  
  
  # Make Heatmap
  png(filename=paste(super_directory, "/", group1,"_vs_",group2,".png",sep=""))
  heatmap.2(matrix_nums_heatmap,
            labRow=matrix_gene_heatmap,
            # more options and specifications
            key=T, trace="none",symkey=T,symbreaks=T,
            Colv=F, Rowv=T,
            # hclustfun=hclust, hclust_method="average",
            col=mypal(20), scale="row", dendrogram ="row",
            cexRow=1, # text size
            margin=c(14, 10), main=title)
  dev.off()
  
  directory <- paste(super_directory, "/", "upregulated_boxplots", sep="")
  dir.create(directory)
  for(gene in up_regulated_20$hgnc_symbol)
  {
    make_boxplot(gene, directory)
  }
  
  directory <- paste(super_directory, "/", "downregulated_boxplots", sep="")
  dir.create(directory)
  for(gene in down_regulated_20$hgnc_symbol)
  {
    make_boxplot(gene, directory)
  }
}
make_boxplot <- function(gene, directory="."){
  plots <- list()
  for(i in unique(group)){
    temp <- as.vector(t(counts[match(gene, as.vector(genes_temp[,2])),group == i]))
    plots[[length(plots)+1]] <- temp
  }
  png(filename=paste(directory, "/", gene, "_Expression.png",sep=""))
  par(mar=c(14, 5, 3, 1))
  boxplot(plots, names=unique(group), las=2)
  mtext(text="Group", side=1, line=12)
  mtext(text="Expression (log Gene Counts)", side=2, line=3)
  mtext(text=paste("Expression of ", gene, "Gene"), side=3, line=1)
  dev.off()
}

# Import data
df <- read.csv("data.csv")
genes_metadata <- read.csv("genes_metadata.csv")
types <- read.csv("types.csv")[,2]
rownames(df) <- make.names(df[,1], unique = TRUE)
all_labels <- read.csv("all_labels.csv")[,2]
df <- df[,c(2:dim(df)[2])]
colnames(df) <- paste(colnames(df), all_labels)
df <- exp(df)
df[c(1:10),c(1:5)] # Take a look at the data

# Palate for Heatmaps
mypal <- colorRampPalette(c("blue", "white", "red"))

# Create DESeq2 object
counts <- df
genes_annotations <- read.csv("annotation.csv")
gene_names <- str_split(genes_metadata$Gene.Symbol, " /// ", simplify=T)[,1]
genes_temp <- cbind(genes_annotations$gene_id[match(tolower(gene_names), tolower(genes_annotations$gene_name))], gene_names)
rownames(genes_temp) <- genes_temp[,1]
rownames(genes_temp)[is.na(rownames(genes_temp))] <- "unknown"
rownames(genes_temp) <- make.names(rownames(genes_temp), unique = TRUE)
rownames(counts) <- rownames(genes_temp)
colnames(genes_temp) <- c("ensembl_gene_id", "hgnc_symbol")
group <- all_labels
all_cds <- DGEList(counts, group=group, genes=genes_temp)
all_keep <- rowSums(edgeR::cpm(all_cds)>3) >= 6
all_count_filtered <- all_cds[all_keep,]
dim(all_cds)[1]
dim(all_count_filtered)[1]

# cpmRNA
cpmRNA <- edgeR::cpm(all_count_filtered)
cpmRNA <- as.data.frame(cpmRNA)

# cf_all, genes_all
cf_all <- calcNormFactors(all_count_filtered, method="TMM")
genes_all <- cf_all$genes$hgnc_symbol
cf_all <- estimateDisp(cf_all) 
cf_all <- estimateCommonDisp(cf_all)
cf_all <- estimateTagwiseDisp(cf_all)

# PCA Visualization
par(bty = 'n')
FCpca <- prcomp(t(edgeR::cpm(cf_all[,types == "OHN"], log=TRUE)))
pcaValues <- data.frame(FCpca$x[,1:2])
pcaValues$group <- group
sds <- FCpca[1]$sdev
variances <- sds^2 / sum(sds^2)
pcaValues
variances

png(filename="pca_OHN.png", width=600, height=600)
group_name <- group[types == "OHN"]
ggplot(pcaValues, aes(PC1, PC2, color = group_name)) +
  geom_point(size = 5) +
  xlab(paste("PC1 (", round(100*variances[1],2), "% Variance)", sep="")) +
  ylab(paste("PC2 (", round(100*variances[2],2), "% Variance)", sep="")) + 
  ggtitle(paste("PCA of RNA-seq OHN (", round(100*(variances[1] + variances[2]),2), "% Variance)", sep=""))
dev.off()

par(bty = 'n')
FCpca <- prcomp(t(edgeR::cpm(cf_all[,types == "retina"], log=TRUE)))
pcaValues <- data.frame(FCpca$x[,1:2])
pcaValues$group <- group
sds <- FCpca[1]$sdev
variances <- sds^2 / sum(sds^2)
pcaValues
variances

png(filename="pca_retina.png", width=600, height=600)
group_name <- group[types == "retina"]
ggplot(pcaValues, aes(PC1, PC2, color = group_name)) +
  geom_point(size = 5) +
  xlab(paste("PC1 (", round(100*variances[1],2), "% Variance)", sep="")) +
  ylab(paste("PC2 (", round(100*variances[2],2), "% Variance)", sep="")) + 
  ggtitle(paste("PCA of RNA-seq retina (", round(100*(variances[1] + variances[2]),2), "% Variance)", sep=""))
dev.off()

# Dispersion Visualization
# dispersion plot
png(filename="dispersion_plot.png")
plotBCV(cf_all)
title("Dispersion Visualization", adj = 0, line = 0.1)
dev.off()
cf_all$common.dispersion
sqrt(cf_all$common.dispersion)

# cpmTMM
cpmTMM<-edgeR::cpm(cf_all)
cpmTMM<-data.frame(cpmTMM)
cpmTMM<-cbind(cpmTMM,cf_all$genes)
cpmTMM<-data.frame(cpmTMM)

# exactTest
unique(all_labels)
freqs <- c()
for(i in unique(all_labels))
{
  freqs <- c(freqs, sum(match(all_labels, i), na.rm=T))
}
cbind(freqs, unique(all_labels))

# Analysis
comparisons <- list(c("Severe OHN", "D2-Gpnmb+ control OHN"), 
                    c("Severe OHN", "Preglaucoma control OHN"),
                    c("Severe OHN", "No or early 1 OHN"),
                    c("Severe retina", "D2-Gpnmb+ control retina"),
                    c("Severe retina", "No or early 1 retina"))
for(i in comparisons)
{
  analysis(i[1], i[2])
}















































































































