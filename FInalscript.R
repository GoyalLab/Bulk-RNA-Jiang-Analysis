# Libraries ####
library(rhdf5) # For hdf5 files from bootstraps.
library(tidyverse) # Provides helpful R packages for data science.
library(tximport) # How Kallisto results are brought into R.
library(biomaRt) #older versions may use ensembldb instead, below
library(EnsDb.Hsapiens.v86) # Human-specific database package.
library(beepr) # Surprise function, very necessary.
library(datapasta) #
library(DESeq2) #differential gene expression analysis
library(ggplot2) #plotting 
library(ggrepel) #labelling plot 
library(pheatmap) # heat maps 
library(dplyr) #good for lots of stuff 
library(org.Hs.eg.db) # human genome data 
library(GOplot) # go plot stuff 
library(scater) #
library(patchwork) #dependency 
library(clusterProfiler) # rnrichment GO analyssi 


# Set Directories ####
dataDirectory <- "/Users/crc2857/Library/CloudStorage/GoogleDrive-cedric.chai123@gmail.com/.shortcut-targets-by-id/1bsXmswJCH-Jf6uILSl2AvBgzP0CbdW7r/JiangRevisions/plotData/"
data2Directory <- "/Users/crc2857/Library/CloudStorage/GoogleDrive-cedric.chai123@gmail.com/.shortcut-targets-by-id/1bsXmswJCH-Jf6uILSl2AvBgzP0CbdW7r/JiangRevisions/extractedData/"
plotDirectory <-  "/Users/crc2857/Library/CloudStorage/GoogleDrive-cedric.chai123@gmail.com/.shortcut-targets-by-id/1bsXmswJCH-Jf6uILSl2AvBgzP0CbdW7r/JiangRevisions/plots/final/final2/"

#Control Data ####
#read in condition data####
diffundiff <- read_tsv(paste0(dataDirectory,"Jiangstudydesign.tsv"))

#bring in file paths for abundance data and confirm exist####
path <- file.path(paste0(data2Directory,diffundiff$...1, "/abundance.tsv"))
all(file.exists(path))

#create Ensembl map in correct form to label data####
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")

#bring in data using file paths and label with gene####
Txi_diffundiff <- tximport(path, type = "kallisto", tx2gene = Tx, txOut = FALSE, countsFromAbundance = "no", ignoreTxVersion = TRUE)

Txi_diffundiffTPM <- tximport(path, type = "kallisto", tx2gene = Tx, txOut = FALSE, countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE)

#round count data to prevent decimals so it can be used for analysis####
Txi_diffundiff$counts= round(Txi_diffundiff$counts)

#separate out just the count and abundance data####
diffundiffcounts <- Txi_diffundiff$counts

diffundiffcounts <- subset(diffundiffcounts, rowSums(diffundiffcounts, na.rm = TRUE) >= 10)

diffundiffcountsTPM <- Txi_diffundiffTPM$abundance

#label columns for count and abundance data####
diffundiffLabels <- diffundiff$...1

colnames(diffundiffcounts) <- c(diffundiffLabels)
colnames(diffundiffcountsTPM) <- c(diffundiffLabels)

#turn condition data into a usable form####
diffundiff <- as.data.frame(diffundiff)

rownames(diffundiff) <- c(diffundiffLabels)

all(rownames(diffundiff) == colnames(diffundiffcounts))

#identify the factor to be analysed####
diffundiff$Differentiated_or_Undifferentiated <- factor(diffundiff$Differentiated_or_Undifferentiated)


#perform DESeq analysis with formatted data####

diffundiffdds <- DESeqDataSetFromMatrix(countData = diffundiffcounts, colData = diffundiff, design = ~ Differentiated_or_Undifferentiated) 

diffundiffdds <- DESeq(diffundiffdds)
diffundiffres <- results(diffundiffdds, contrast=c("Differentiated_or_Undifferentiated", "Differentiated", "Undifferentiated"))

#create PCA plot####
diffundiffvsd <- vst(diffundiffdds, blind=FALSE)

plotdiffundiff <- plotPCA(diffundiffvsd, intgroup=c("Differentiated_or_Undifferentiated"), returnData=TRUE)

#percentVar <- round(100 * attr(plotdiffundiff, "percentVar"))
#ggplot(plotdiffundiff) + aes(PC1, PC2, colour=Differentiated_or_Undifferentiated) + geom_text_repel(label = rownames(plotdiffundiff)) + geom_point(size=3) +theme_bw() + labs(title= "Differentiated Vs Unifferentiated", subtitle = "Plot1") + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed() 

#plot1
percentVar <- round(100 * attr(plotdiffundiff, "percentVar"))
plotcontrol = ggplot(plotdiffundiff) + 
  aes(PC1, PC2, colour=Differentiated_or_Undifferentiated) + 
  geom_text_repel(label = rownames(plotdiffundiff)) + 
  geom_point(size=3) +
  theme_classic() + 
  labs(title= "Differentiated Vs Undifferentiated", subtitle = "Plotcontrol") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), line = element_blank(), legend.position= "none") +
  coord_fixed() 

ggsave(plotcontrol, file = paste0(plotDirectory, 'PCA_Control_Wtext.svg'), width = 6, height = 4)

plotcontrola = ggplot(plotdiffundiff) + 
  aes(PC1, PC2, colour=Differentiated_or_Undifferentiated) + 
  geom_point(size=3) +
  theme_classic() + 
  labs(title= "Differentiated Vs Undifferentiated", subtitle = "Plotcontrol") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), line = element_blank(), legend.position= "none") +
  coord_fixed() 
ggsave(plotcontrola, file = paste0(plotDirectory, 'PCA_Control.svg'), width = 6, height = 4)


#**
diffundiffSig <- subset(diffundiffres, padj < 0.01)
diffundiffres.df <- as.data.frame(diffundiffSig)

#**
topdiffundiffres.df <- diffundiffres.df[(diffundiffres.df$baseMean > 10) & (abs(diffundiffres.df$log2FoldChange) > 2), ]
topdiffundiffres.df %>% drop_na(log2FoldChange)

topdiffundiffres.df <- topdiffundiffres.df[order(topdiffundiffres.df$log2FoldChange, decreasing = TRUE),]

topdiffundiffres.df$gene <- rownames(topdiffundiffres.df)

num_keep <- 50

topdiffundiffrowkeep <- c(seq(1:num_keep), seq((nrow(topdiffundiffres.df)-num_keep), nrow(topdiffundiffres.df)))

alldiffundiff <- merge(diffundiffcountsTPM, topdiffundiffres.df, by = 0)
alldiffundiffcounts <- alldiffundiff[,2:21]
rownames(alldiffundiffcounts) <- topdiffundiffres.df$gene

filteredalldiffundiffcounts <- as.data.frame(alldiffundiffcounts[topdiffundiffrowkeep,])

controlheatmap = pheatmap(filteredalldiffundiffcounts, scale = "row", main = "All Samples", fontsize_row = 6)
ggsave(controlheatmap, file = paste0(plotDirectory, 'Heatmap_Control_Wtext.svg'), width = 5, height = 10)

mycolours <- c("blue", "red", "black")
names(mycolours) <- c("Down regulated", "Up regulated", "Nonsignificant")

subdiffundiffresOrdered <- as.data.frame(diffundiffres[order(diffundiffres$pvalue),])
subdiffundiffresOrdered %>% drop_na(log2FoldChange)
subdiffundiffresOrdered %>% drop_na(padj)

subdiffundiffresOrdered$differentialexpression <- "Nonsignificant"
subdiffundiffresOrdered$differentialexpression[subdiffundiffresOrdered$log2FoldChange > 2 & subdiffundiffresOrdered$padj < 0.01] <- "Up regulated"
subdiffundiffresOrdered$differentialexpression[subdiffundiffresOrdered$log2FoldChange < -2 & subdiffundiffresOrdered$padj < 0.01] <- "Down regulated"
subdiffundiffresOrdered$deLabel <- NA
subdiffundiffresOrdered$deLabel[subdiffundiffresOrdered$differentialexpression != "Nonsignificant"] <- rownames(subdiffundiffresOrdered)[subdiffundiffresOrdered$differentialexpression != "Nonsignificant"]
#volcano plot
controlvolcano =  ggplot(subdiffundiffresOrdered, aes(x= log2FoldChange, y= -log10(padj), col = differentialexpression, label = deLabel)) + geom_point() + geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed") + geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") + geom_text_repel() + scale_color_manual(values = mycolours) + labs(title= " Differentiated Vs Undifferentiated") + xlab(bquote(log[2]*FoldChange)) + ylab(bquote(-log[10]*(padj))) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), line = element_blank(), legend.position= "none") +
  coord_fixed() 
ggsave(controlvolcano, file = paste0(plotDirectory, 'Volcano_Control_Wtext.svg'), width = 6, height = 4)

controlvolcanoa =  ggplot(subdiffundiffresOrdered, aes(x= log2FoldChange, y= -log10(padj), col = differentialexpression, label = deLabel)) + geom_point() + geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed") + geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") + geom_text_repel() + scale_color_manual(values = mycolours) + labs(title= " Differentiated Vs Undifferentiated") + xlab(bquote(log[2]*FoldChange)) + ylab(bquote(-log[10]*(padj))) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), line = element_blank(), legend.position= "none") +
  coord_fixed() 
ggsave(controlvolcano, file = paste0(plotDirectory, 'Volcano_Control.svg'), width = 6, height = 4)



#**
diffundiffresNONA <- as.data.frame(diffundiffres)
diffundiffresNONA %>% drop_na(log2FoldChange)
diffundiffresNONA %>% drop_na(padj)
updiffundiffDGE <- diffundiffresNONA[(diffundiffresNONA$log2FoldChange > 2) & (diffundiffresNONA$padj < 0.05), ]
downdiffundiffDGE <- diffundiffresNONA[(diffundiffresNONA$log2FoldChange < -2) & (diffundiffresNONA$padj < 0.05), ]


updiffundiffDGE$gene_name <- rownames(updiffundiffDGE)
downdiffundiffDGE$gene_name <- rownames(downdiffundiffDGE)

myMart <- useMart(biomart= "ENSEMBL_MART_ENSEMBL")
available.packages <- listDatasets(myMart)
human.anno <- useMart(biomart= "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
human.attributes <- listFilters(human.anno)
humanentrezgene <- getBM( attributes = c("ensembl_transcript_id_version", "external_gene_name", "entrezgene_id"), mart = human.anno, values = "with_entrezgene")
humanentrezgene <- as.data.frame(humanentrezgene)
humanentrezgene <- dplyr::rename(humanentrezgene, target_id = ensembl_transcript_id_version, gene_name = external_gene_name )

updiffundiffDGEjoin <- merge( x = updiffundiffDGE, y = humanentrezgene, by = "gene_name", all.x =TRUE )
updiffundiffDGEjoin <- updiffundiffDGEjoin %>% drop_na()

#**
updiffundiffGO <- enrichGO( gene = updiffundiffDGEjoin$entrezgene_id, OrgDb = org.Hs.eg.db, ont = "all", qvalueCutoff = 0.01)

downdiffundiffDGEjoin <- merge( x = downdiffundiffDGE, y = humanentrezgene, by = "gene_name", all.x =TRUE )
downdiffundiffDGEjoin <- downdiffundiffDGEjoin %>% drop_na()

#**
downdiffundiffGO <- enrichGO( gene = downdiffundiffDGEjoin$entrezgene_id, OrgDb = org.Hs.eg.db, ont = "all", qvalueCutoff = 0.01)

upbardiffundiffGO =  barplot(updiffundiffGO, showCategory = 15) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=15), 
        legend.key.size = unit(0.6, 'cm'), legend.title = element_text(size=12), legend.text = element_text(size=12)) + 
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(upbardiffundiffGO, file = paste0(plotDirectory, 'Up_GO_Control_Wtext.svg'), width = 8, height = 7)

upbardiffundiffGOa =  barplot(updiffundiffGO, showCategory = 15) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=15), 
        legend.key.size = unit(0.6, 'cm'), legend.title = element_text(size=12), legend.text = element_text(size=12)) + 
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(upbardiffundiffGOa, file = paste0(plotDirectory, 'Up_GO_Control.svg'), width = 8, height = 7)



downbardiffundiffGO = barplot(downdiffundiffGO, showCategory = 15) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=15), 
        legend.key.size = unit(0.6, 'cm'), legend.title = element_text(size=12), legend.text = element_text(size=12)) + 
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(downbardiffundiffGO, file = paste0(plotDirectory, 'Down_GO_Control_Wtext.svg'), width = 8, height = 7)


downbardiffundiffGOa = barplot(downdiffundiffGO, showCategory = 15) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=15), 
        legend.key.size = unit(0.6, 'cm'), legend.title = element_text(size=12), legend.text = element_text(size=12)) + 
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggsave(downbardiffundiffGOa, file = paste0(plotDirectory, 'Down_GO_Control.svg'), width = 8, height = 7)

nrow(topdiffundiffres.df)
#3464
#Grouping Data ####
#read in condition data####
healthy <- read_tsv(paste0(dataDirectory,"Healthyonly.tsv"))
diseased <- read_tsv(paste0(dataDirectory,"Diseasedonly.tsv"))
differentiated <- read_tsv(paste0(dataDirectory,"Differentiatedonly.tsv"))
undifferentiated <- read_tsv(paste0(dataDirectory,"Undifferentiatedonly.tsv"))

#bring in file paths for abundance data and confirm exist####
differentiatedpath <- file.path(paste0(data2Directory,differentiated$...1, "/abundance.tsv"))
healthypath <- file.path(paste0(data2Directory,healthy$...1, "/abundance.tsv"))
undifferentiatedpath <- file.path(paste0(data2Directory,undifferentiated$...1, "/abundance.tsv"))
diseasedpath <- file.path(paste0(data2Directory,diseased$...1, "/abundance.tsv"))
all(file.exists(differentiatedpath))
all(file.exists(undifferentiatedpath))
all(file.exists(diseasedpath))
all(file.exists(healthypath))


#bring in data using file paths and label with gene####
Txi_differentiated <- tximport(differentiatedpath, type = "kallisto", tx2gene = Tx, txOut = FALSE, countsFromAbundance = "no", ignoreTxVersion = TRUE)
Txi_undifferentiated <- tximport(undifferentiatedpath, type = "kallisto", tx2gene = Tx, txOut = FALSE, countsFromAbundance = "no", ignoreTxVersion = TRUE)
Txi_diseased <- tximport(diseasedpath, type = "kallisto", tx2gene = Tx, txOut = FALSE, countsFromAbundance = "no", ignoreTxVersion = TRUE)
Txi_healthy <- tximport(healthypath, type = "kallisto", tx2gene = Tx, txOut = FALSE, countsFromAbundance = "no", ignoreTxVersion = TRUE)

Txi_differentiatedTPM <- tximport(differentiatedpath, type = "kallisto", tx2gene = Tx, txOut = FALSE, countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE)
Txi_undifferentiatedTPM <- tximport(undifferentiatedpath, type = "kallisto", tx2gene = Tx, txOut = FALSE, countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE)
Txi_diseasedTPM <- tximport(diseasedpath, type = "kallisto", tx2gene = Tx, txOut = FALSE, countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE)
Txi_healthyTPM <- tximport(healthypath, type = "kallisto", tx2gene = Tx, txOut = FALSE, countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE)

#Round counts to integers.
Txi_differentiated$counts= round(Txi_differentiated$counts)
Txi_undifferentiated$counts= round(Txi_undifferentiated$counts)
Txi_diseased$counts= round(Txi_diseased$counts)
Txi_healthy$counts= round(Txi_healthy$counts)

differentiatedcounts <- Txi_differentiated$counts #55957
undifferentiatedcounts <- Txi_undifferentiated$counts
diseasedcounts <- Txi_diseased$counts
healthycounts <- Txi_healthy$counts

differentiatedcountsTPM <- Txi_differentiatedTPM$abundance #55957
undifferentiatedcountsTPM <- Txi_undifferentiatedTPM$abundance
diseasedcountsTPM <- Txi_diseasedTPM$abundance
healthycountsTPM <- Txi_healthyTPM$abundance

#label columns for count data####
differentiatedLabels <- differentiated$...1
undifferentiatedLabels <- undifferentiated$...1
diseasedLabels <- diseased$...1
healthyLabels <- healthy$...1

colnames(undifferentiatedcounts) <- c(undifferentiatedLabels)
colnames(diseasedcounts) <- c(diseasedLabels)
colnames(healthycounts) <- c(healthyLabels)
colnames(differentiatedcounts) <- c(differentiatedLabels)

differentiatedcounts <- subset(differentiatedcounts, rowSums(differentiatedcounts, na.rm = TRUE) >= 10)
undifferentiatedcounts <- subset(undifferentiatedcounts, rowSums(undifferentiatedcounts, na.rm = TRUE) >= 10) 
diseasedcounts <- subset(diseasedcounts, rowSums(diseasedcounts, na.rm = TRUE) >= 10) 
healthycounts <- subset(healthycounts, rowSums(healthycounts, na.rm = TRUE) >= 10) 

colnames(undifferentiatedcountsTPM) <- c(undifferentiatedLabels)
colnames(diseasedcountsTPM) <- c(diseasedLabels)
colnames(healthycountsTPM) <- c(healthyLabels)
colnames(differentiatedcountsTPM) <- c(differentiatedLabels)

#turn condition data into a usable form####
healthy <- as.data.frame(healthy)
diseased <- as.data.frame((diseased))
differentiated <- as.data.frame(differentiated)
undifferentiated <- as.data.frame(undifferentiated)

rownames(healthy) <- c(healthyLabels)
rownames(diseased) <- c(diseasedLabels)
rownames(differentiated) <- c(differentiatedLabels)
rownames(undifferentiated) <- c(undifferentiatedLabels)

all(rownames(differentiated) == colnames(differentiatedcounts))
all(rownames(healthy) == colnames(healthycounts))
all(rownames(diseased)  == colnames(diseasedcounts))
all(rownames(undifferentiated) == colnames(undifferentiatedcounts))

#identify the factor to be analysed####
differentiated$Healthy_or_Disease <- factor(differentiated$Healthy_or_Disease)
undifferentiated$Healthy_or_Disease <- factor(undifferentiated$Healthy_or_Disease)
diseased$Differentiated_or_Undifferentiated <- factor(diseased$Differentiated_or_Undifferentiated)
healthy$Differentiated_or_Undifferentiated <- factor(healthy$Differentiated_or_Undifferentiated)

#perform DESeq analysis with formatted data####
differentiateddds <- DESeqDataSetFromMatrix(countData = differentiatedcounts, colData = differentiated, design = ~ Healthy_or_Disease)
undifferentiateddds <- DESeqDataSetFromMatrix(countData = undifferentiatedcounts, colData = undifferentiated, design = ~ Healthy_or_Disease)
healthyddds <- DESeqDataSetFromMatrix(countData = healthycounts, colData = healthy, design = ~ Differentiated_or_Undifferentiated)
diseasedddds <- DESeqDataSetFromMatrix(countData = diseasedcounts, colData = diseased, design = ~ Differentiated_or_Undifferentiated) 

diseasedddds <- DESeq(diseasedddds)
diseasedres <- results(diseasedddds)

healthyddds <- DESeq(healthyddds)
healthyres <- results(healthyddds)

differentiateddds <- DESeq(differentiateddds)
differentiatedres <- results(differentiateddds)

undifferentiateddds <- DESeq(undifferentiateddds)
undifferentiatedres <- results(undifferentiateddds)

#create PCA plot####
differentiatedvsd <- vst(differentiateddds, blind=FALSE)
undifferentiatedvsd <- vst(undifferentiateddds, blind = FALSE)
diseasedvsd <- vst(diseasedddds, blind = FALSE)
healthyvsd <- vst(healthyddds, blind = FALSE)

plothealthy <- plotPCA(healthyvsd, intgroup=c("Differentiated_or_Undifferentiated"), returnData=TRUE)
plotdiseased <- plotPCA(diseasedvsd, intgroup=c("Differentiated_or_Undifferentiated"), returnData=TRUE)
plotundifferentiated <- plotPCA(undifferentiatedvsd, intgroup=c("Healthy_or_Disease"), returnData=TRUE) 
plotdifferentiated <- plotPCA(differentiatedvsd, intgroup=c("Healthy_or_Disease"), returnData=TRUE)

#plot1
percentVar <- round(100 * attr(plotdifferentiated, "percentVar"))
plot1 = ggplot(plotdifferentiated) + 
  aes(PC1, PC2, colour=Healthy_or_Disease) + 
  geom_text_repel(label = rownames(plotdifferentiated)) + 
  geom_point(size=3) +
  theme_classic() + 
  labs(title= "Healthy Differentiated Vs Diseased Differentiated", subtitle = "Plot1") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), line = element_blank(), legend.position= "none") +
  coord_fixed() 

ggsave(plot1, file = paste0(plotDirectory, 'PCA_Differentiated_Wtext.svg'), width = 6, height = 4)

plot1a = ggplot(plotdifferentiated) + 
  aes(PC1, PC2, colour=Healthy_or_Disease) + 
  geom_point(size=3) +
  theme_classic() + 
  labs(title= "Healthy Differentiated Vs Diseased Differentiated", subtitle = "Plot1") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), line = element_blank(), legend.position= "none") +
  coord_fixed() 
ggsave(plot1a, file = paste0(plotDirectory, 'PCA_Differentiated.svg'), width = 6, height = 4)

#plot2
percentVar <- round(100 * attr(plothealthy, "percentVar"))

plot2 = ggplot(plothealthy) + 
  aes(PC1, PC2, colour=Differentiated_or_Undifferentiated) + 
  geom_text_repel(label = rownames(plothealthy)) + 
  geom_point(size=3) + 
  theme_classic() + 
  labs(title= "Healthy Differentiated Vs Healthy Undifferentiated", subtitle = "Plot2") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +     
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +  
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), line = element_blank(), legend.position= "none") +
  coord_fixed() 
ggsave(plot2, file = paste0(plotDirectory, 'PCA_Healthy_Wtext.svg'), width = 6, height = 4)

plot2a = ggplot(plothealthy) + 
  aes(PC1, PC2, colour=Differentiated_or_Undifferentiated) + 
  geom_point(size=3) + 
  theme_classic() + 
  labs(title= "Healthy Differentiated Vs Healthy Undifferentiated", subtitle = "Plot2") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +     
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +  
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), line = element_blank(), legend.position= "none") +
  coord_fixed() 
ggsave(plot2a, file = paste0(plotDirectory, 'PCA_Healthy.svg'), width = 6, height = 4)

#plot3
percentVar <- round(100 * attr(plotundifferentiated, "percentVar"))
plot3 = ggplot(plotundifferentiated, aes(PC1, PC2, color=Healthy_or_Disease)) +
  geom_point(size=3) +
  theme_classic() + 
  geom_text_repel(label = rownames(plotundifferentiated)) + 
  labs(title= "Healthy Undifferentiated Vs Diseased Undifferentiated ", subtitle = "Plot3") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +  
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), line = element_blank(), legend.position= "none") +
  coord_fixed() 
ggsave(plot3, file = paste0(plotDirectory, 'PCA_undifferentiated_Wtext.svg'), width = 6, height = 4)

plot3a = ggplot(plotundifferentiated, aes(PC1, PC2, color=Healthy_or_Disease)) +
  geom_point(size=3) +
  theme_classic() + 
  labs(title= "Healthy Undifferentiated Vs Diseased Undifferentiated ", subtitle = "Plot3") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +  
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), line = element_blank(), legend.position= "none") +
  coord_fixed() 
ggsave(plot3a, file = paste0(plotDirectory, 'PCA_undifferentiated.svg'), width = 6, height = 4)

#plot4
percentVar <- round(100 * attr(plotdiseased, "percentVar"))
plot4 = ggplot(plotdiseased, aes(PC1, PC2, color=Differentiated_or_Undifferentiated)) + 
  geom_point(size=3) + 
  theme_classic() + 
  geom_text_repel(label = rownames(plotdiseased)) + 
  labs(title= "Diseased Differentiated Vs Diseased Undifferentiated", subtitle = "Plot4") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), line = element_blank(), legend.position= "none") +
  coord_fixed()

ggsave(plot4, file = paste0(plotDirectory, 'PCA_disease_Wtext.svg'), width = 6, height = 4)

plot4a = ggplot(plotdiseased, aes(PC1, PC2, color=Differentiated_or_Undifferentiated)) + 
  geom_point(size=3) + 
  theme_classic() + 
  geom_text_repel(label = rownames(plotdiseased)) + 
  labs(title= "Diseased Differentiated Vs Diseased Undifferentiated", subtitle = "Plot4") + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), line = element_blank(), legend.position= "none") +
  coord_fixed()

ggsave(plot4a, file = paste0(plotDirectory, 'PCA_disease.svg'), width = 6, height = 4)


#first filter ####
# **Note thresholds of p adjusted value less than 0.01
healthyresSig <- subset(healthyres, padj < 0.01)
diseasedresSig <- subset(diseasedres, padj < 0.01)
differentiatedresSig <- subset(differentiatedres, padj < 0.01)
undifferentiatedresSig <- subset(undifferentiatedres, padj < 0.01)

#turn into data frame
undifferentiatedres.df <- as.data.frame(undifferentiatedresSig)
differentiatedres.df <- as.data.frame(differentiatedresSig)
healthyres.df <- as.data.frame(healthyresSig)
diseasedres.df <- as.data.frame(diseasedresSig)

#second filter, removing insignificant data and genes with NA values for log2 fold change
topundifferentiatedres.df <- undifferentiatedres.df[(undifferentiatedres.df$baseMean > 10) & (abs(undifferentiatedres.df$log2FoldChange) > 2),] 
topundifferentiatedres.df %>% drop_na(log2FoldChange)

topdifferentiatedres.df <- differentiatedres.df[(differentiatedres.df$baseMean > 10) & (abs(differentiatedres.df$log2FoldChange) > 2), ]
topdifferentiatedres.df %>% drop_na(log2FoldChange)

tophealthyres.df <- healthyres.df[(healthyres.df$baseMean > 10) & (abs(healthyres.df$log2FoldChange) > 2),]
tophealthyres.df %>% drop_na(log2FoldChange)

topdiseasedres.df <- diseasedres.df[(diseasedres.df$baseMean) > 10 & (abs(diseasedres.df$log2FoldChange) > 2),]
topdiseasedres.df %>% drop_na(log2FoldChange)

#ordering the data to have the top log2fold genes on top and the lowest non NA values at the bottom
topundifferentiatedres.df <- topundifferentiatedres.df[order(topundifferentiatedres.df$log2FoldChange, decreasing = TRUE),]
topdifferentiatedres.df <- topdifferentiatedres.df[order(topdifferentiatedres.df$log2FoldChange, decreasing = TRUE),]
tophealthyres.df <- tophealthyres.df[order(tophealthyres.df$log2FoldChange, decreasing = TRUE),]
topdiseasedres.df <- topdiseasedres.df[order(topdiseasedres.df$log2FoldChange, decreasing = TRUE),]

#create a gene row so its easier to search
topundifferentiatedres.df$gene <- rownames(topundifferentiatedres.df)
topdifferentiatedres.df$gene <- rownames(topdifferentiatedres.df)
topdiseasedres.df$gene <- rownames(topdiseasedres.df)
tophealthyres.df$gene <- rownames(tophealthyres.df)

#DGE list 
write.csv(as.data.frame(tophealthyres.df), 
         file= paste0(plotDirectory, "Healthydifferentiated_Healthyundifferentiated.csv"))

write.csv(as.data.frame(topdiseasedres.df), 
         file= paste0(plotDirectory,"Diseaseddifferentiated_Diseasedundifferentiated.csv"))

write.csv(as.data.frame(topdifferentiatedres.df), 
          file= paste0(plotDirectory,"healthyDifferentiated_diseasedDifferentiated.csv"))

write.csv(as.data.frame(topundifferentiatedres.df), 
          file= paste0(plotDirectory,"healthyUnifferentiated_diseasedUnifferentiated.csv"))

write.csv(as.data.frame(topdiffundiffres.df), 
          file= paste0(plotDirectory, "Control.csv"))


# since the data is organised by log numbers we can find the row numbers of the 50 most up and down regulated genes
num_keep <- 50
topdifferentiatedrowkeep <- c(seq(1:num_keep), seq((nrow(topdifferentiatedres.df)-num_keep), nrow(topdifferentiatedres.df)))
topundifferentiatedrowkeep <- c(seq(1:num_keep), seq((nrow(topundifferentiatedres.df)-num_keep), nrow(topundifferentiatedres.df)))
topdiseasedrowkeep <- c(seq(1:num_keep), seq((nrow(topdiseasedres.df)-num_keep), nrow(topdiseasedres.df)))
tophealthyrowkeep <- c(seq(1:num_keep), seq((nrow(tophealthyres.df)-num_keep), nrow(tophealthyres.df)))


topundifferentiatedres.df$gene <- rownames(topundifferentiatedres.df)
topdifferentiatedres.df$gene <- rownames(topdifferentiatedres.df)
topdiseasedres.df$gene <- rownames(topdiseasedres.df)
tophealthyres.df$gene <- rownames(tophealthyres.df)

# merging the count data will filter out some insignificant  genes
allundifferentiated <- merge(undifferentiatedcountsTPM, topundifferentiatedres.df, by = 0)
alldifferentiated <- merge(differentiatedcountsTPM, topdifferentiatedres.df, by = 0)
alldiseased <- merge(diseasedcountsTPM, topdiseasedres.df, by = 0)
allhealthy <- merge(healthycountsTPM, tophealthyres.df, by = 0)

#the log data was useful for filtering but the only information to make the heat map is count data, sample name and gene name. Remove anything else.
allundifferentiatedcountsTPM <- allundifferentiated[,2:11]
alldifferentiatedcountsTPM <- alldifferentiated[,2:11]
allhealthycountsTPM <- allhealthy[,2:11]
alldiseasedcountsTPM <- alldiseased[,2:11]

rownames(allundifferentiatedcountsTPM) <- topundifferentiatedres.df$gene
rownames(alldifferentiatedcountsTPM) <- topdifferentiatedres.df$gene
rownames(alldiseasedcountsTPM) <- topdiseasedres.df$gene
rownames(allhealthycountsTPM) <- tophealthyres.df$gene

filteredallhealthycounts <- as.data.frame(allhealthycountsTPM[tophealthyrowkeep,])
filteredalldiseasedcounts <- as.data.frame(alldiseasedcountsTPM[topdiseasedrowkeep,])
filteredalldifferentiatedcounts <- as.data.frame(alldifferentiatedcountsTPM[topdifferentiatedrowkeep,])
filteredallundifferentiatedcounts <- as.data.frame(allundifferentiatedcountsTPM)

library(pheatmap)
#heatmap1
heatmap1 <- pheatmap(filteredallundifferentiatedcounts, scale = "row", main = "Undifferentiated Samples", fontsize_row = 6)
ggsave(heatmap1, file = paste0(plotDirectory, 'Heatmap1.svg'), width = 5, height = 10)

#heatmap2
heatmap2 <- pheatmap(filteredalldifferentiatedcounts, scale = "row", main = "Differentiated Samples", fontsize_row = 6)
ggsave(heatmap2, file = paste0(plotDirectory, 'Heatmap2.svg'), width = 5, height = 10)

#heatmap3
heatmap3 <- pheatmap(filteredalldiseasedcounts, scale = "row", main = "Diseased Samples", fontsize_row = 6)
ggsave(heatmap3, file = paste0(plotDirectory, 'Heatmap3.svg'), width = 5, height = 10)

#heatmap4
heatmap4 <- pheatmap(filteredallhealthycounts, scale = "row", main = "Healthy Samples", fontsize_row = 6)
ggsave(heatmap4, file = paste0(plotDirectory, 'Heatmap4.svg'), width = 5, height = 10)

subhealthyresOrdered <- as.data.frame(healthyres[order(healthyres$pvalue),])
subdiseasedresOrdered <- as.data.frame(diseasedres[order(diseasedres$pvalue),])
subdifferentiatedresOrdered <- as.data.frame(differentiatedres[order(differentiatedres$pvalue),])
subundifferentiatedresOrdered <- as.data.frame(undifferentiatedres[order(undifferentiatedres$pvalue),])

#remove genes that have NA values for log2fold change or adjusted p value
subhealthyresOrdered %>% drop_na(log2FoldChange)
subhealthyresOrdered %>% drop_na(padj)

subdiseasedresOrdered %>% drop_na(log2FoldChange)
subdiseasedresOrdered %>% drop_na(padj)

subdifferentiatedresOrdered %>% drop_na(log2FoldChange)
subdifferentiatedresOrdered %>% drop_na(padj)

subundifferentiatedresOrdered %>% drop_na(log2FoldChange)
subundifferentiatedresOrdered %>% drop_na(padj)

#creating a volcano plot with labels for up and down regulated genes####
subhealthyresOrdered$differentialexpression <- "Nonsignificant"
subhealthyresOrdered$differentialexpression[subhealthyresOrdered$log2FoldChange > 2 & subhealthyresOrdered$padj < 0.05] <- "Up regulated"
subhealthyresOrdered$differentialexpression[subhealthyresOrdered$log2FoldChange < -2 & subhealthyresOrdered$padj < 0.05] <- "Down regulated"
subhealthyresOrdered$deLabel <- NA
subhealthyresOrdered$deLabel[subhealthyresOrdered$differentialexpression != "Nonsignificant"] <- rownames(subhealthyresOrdered)[subhealthyresOrdered$differentialexpression != "Nonsignificant"]
#volcano1
volcano1 <- ggplot(subhealthyresOrdered, aes(x= log2FoldChange, y= -log10(padj), col = differentialexpression, label = deLabel)) + geom_point() + geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed") + geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") + geom_text_repel() + scale_color_manual(values = mycolours) + labs(title= "Healthy Differentiated Vs Healthy Undifferentiated") + xlab(bquote(log[2]*FoldChange)) + ylab(bquote(-log[10]*(padj))) +
theme(panel.border = element_rect(colour = "black", fill=NA, size=1), line = element_blank(), legend.position= "none") +
coord_fixed() 
ggsave(volcano1, file = paste0(plotDirectory, 'Volcano_Healthy.svg'), width = 6, height = 4)

subdiseasedresOrdered$differentialexpression <- "Nonsignificant"
subdiseasedresOrdered$differentialexpression[subdiseasedresOrdered$log2FoldChange > 2 & subdiseasedresOrdered$padj < 0.05] <- "Up regulated"
subdiseasedresOrdered$differentialexpression[subdiseasedresOrdered$log2FoldChange < -2 & subdiseasedresOrdered$padj < 0.05] <- "Down regulated"
subdiseasedresOrdered$deLabel <- NA
subdiseasedresOrdered$deLabel[subdiseasedresOrdered$differentialexpression != "Nonsignificant"] <- rownames(subdiseasedresOrdered)[subdiseasedresOrdered$differentialexpression != "Nonsignificant"]
#volcano2
volcano2 <- ggplot(subdiseasedresOrdered, aes(x= log2FoldChange, y= -log10(padj), col = differentialexpression, label = deLabel)) + geom_point() + geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed") + geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") + geom_text_repel() + scale_color_manual(values = mycolours)+ labs(title= "Healthy Undifferentiated Vs Diseased Undifferentiated ") + xlab(bquote(log[2]*FoldChange)) + ylab(bquote(-log[10]*(padj))) +
theme(panel.border = element_rect(colour = "black", fill=NA, size=1), line = element_blank(), legend.position= "none") +
coord_fixed() 
ggsave(volcano2, file = paste0(plotDirectory, 'Volcano_Diseased.svg'), width = 6, height = 4)

subdifferentiatedresOrdered$differentialexpression <- "Nonsignificant"
subdifferentiatedresOrdered$differentialexpression[subdifferentiatedresOrdered$log2FoldChange > 2 & subdifferentiatedresOrdered$padj < 0.05] <- "Up regulated"
subdifferentiatedresOrdered$differentialexpression[subdifferentiatedresOrdered$log2FoldChange < -2 & subdifferentiatedresOrdered$padj < 0.05] <- "Down regulated"
subdifferentiatedresOrdered$deLabel <- NA
subdifferentiatedresOrdered$deLabel[subdifferentiatedresOrdered$differentialexpression != "Nonsignificant"] <- rownames(subdifferentiatedresOrdered)[subdifferentiatedresOrdered$differentialexpression != "Nonsignificant"]
#volcano3
volcano3 <- ggplot(subdifferentiatedresOrdered, aes(x= log2FoldChange, y= -log10(padj), col = differentialexpression, label = deLabel)) + geom_point() + geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed") + geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") + geom_text_repel() + scale_color_manual(values = mycolours) + labs(title= "Diseased Differentiated Vs Diseased Undifferentiated") + xlab(bquote(log[2]*FoldChange)) + ylab(bquote(-log[10]*(padj))) +
theme(panel.border = element_rect(colour = "black", fill=NA, size=1), line = element_blank(), legend.position= "none") +
coord_fixed() 
ggsave(volcano3, file = paste0(plotDirectory, 'Volcano_Differentiated.svg'), width = 6, height = 4)

subundifferentiatedresOrdered$differentialexpression <- "Nonsignificant"
subundifferentiatedresOrdered$differentialexpression[subundifferentiatedresOrdered$log2FoldChange > 2 & subundifferentiatedresOrdered$padj < 0.05] <- "Up regulated"
subundifferentiatedresOrdered$differentialexpression[subundifferentiatedresOrdered$log2FoldChange < -2 & subundifferentiatedresOrdered$padj < 0.05] <- "Down regulated"
subundifferentiatedresOrdered$deLabel <- NA
subundifferentiatedresOrdered$deLabel[subundifferentiatedresOrdered$differentialexpression != "Nonsignificant"] <- rownames(subundifferentiatedresOrdered)[subundifferentiatedresOrdered$differentialexpression != "Nonsignificant"]
#volcano4
volcano4 <- ggplot(subundifferentiatedresOrdered, aes(x= log2FoldChange, y= -log10(padj), col = differentialexpression, label = deLabel)) + geom_point() + geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed") + geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") + geom_text_repel() + scale_color_manual(values = mycolours) + labs(title= "Healthy Undifferentiated Vs Diseased Undifferentiated ") + xlab(bquote(log[2]*FoldChange)) + ylab(bquote(-log[10]*(padj))) +
theme(panel.border = element_rect(colour = "black", fill=NA, size=1), line = element_blank(), legend.position= "none") +
coord_fixed() 
ggsave(volcano4, file = paste0(plotDirectory, 'Volcano_Undifferentiated.svg'), width = 6, height = 4)

#Make DGE filtered List for Go analysis####
undifferentiatedresNONA <- as.data.frame(undifferentiatedres)
undifferentiatedresNONA %>% drop_na(log2FoldChange)
undifferentiatedresNONA %>% drop_na(padj)
undifferentiatedDGE <- undifferentiatedresNONA[(abs(undifferentiatedresNONA$log2FoldChange) > 2) & (undifferentiatedresNONA$padj < 0.05),] 
undifferentiatedDGE$gene_name <- rownames(undifferentiatedDGE)

differentiatedresNONA <- as.data.frame(differentiatedres)
differentiatedresNONA %>% drop_na(log2FoldChange)
differentiatedresNONA %>% drop_na(padj)
differentiatedDGE <- differentiatedresNONA[(abs(differentiatedresNONA$log2FoldChange) > 2) & (differentiatedresNONA$padj < 0.05), ]
differentiatedDGE$gene_name <- rownames(differentiatedDGE)

healthyresNONA <- as.data.frame(healthyres)
healthyresNONA %>% drop_na(log2FoldChange)
healthyresNONA %>% drop_na(padj)
healthyDGE <- healthyresNONA[(abs(healthyresNONA$log2FoldChange) > 2) & (healthyresNONA$padj < 0.05),]
healthyDGE$gene_name <- rownames(healthyDGE)

diseasedresNONA <- as.data.frame(diseasedres)
diseasedresNONA %>% drop_na(log2FoldChange)
diseasedresNONA %>% drop_na(padj)
diseasedDGE <- diseasedresNONA[(abs(diseasedresNONA$log2FoldChange) > 2) & (diseasedresNONA$padj < 0.05),]
diseasedDGE$gene_name <- rownames(diseasedDGE)

#merging files creates a filtered data set with the entrezgene ID and DESeq data.
undifferentiatedDGEjoin <- merge( x = undifferentiatedDGE, y = humanentrezgene, by = "gene_name", all.x =TRUE )
undifferentiatedDGEjoin <- undifferentiatedDGEjoin %>% drop_na()

differentiatedDGEjoin <- merge( x = differentiatedDGE, y = humanentrezgene, by = "gene_name", all.x =TRUE )
differentiatedDGEjoin <- differentiatedDGEjoin %>% drop_na()

healthyDGEjoin <- merge( x = healthyDGE, y = humanentrezgene, by = "gene_name", all.x =TRUE )
healthyDGEjoin <- healthyDGEjoin %>% drop_na()

diseasedDGEjoin <- merge( x = diseasedDGE, y = humanentrezgene, by = "gene_name", all.x =TRUE )
diseasedDGEjoin <- diseasedDGEjoin %>% drop_na()

# Running GO analysis####
undifferentiatedGO <- enrichGO( gene = undifferentiatedDGEjoin$entrezgene_id, OrgDb = org.Hs.eg.db, ont = "all", qvalueCutoff = 0.01)

differentiatedGO <- enrichGO( gene = differentiatedDGEjoin$entrezgene_id, OrgDb = org.Hs.eg.db, ont = "all", qvalueCutoff = 0.01)

healthyGO <- enrichGO( gene = healthyDGEjoin$entrezgene_id, OrgDb = org.Hs.eg.db, ont = "all", qvalueCutoff = 0.01)

diseasedGO <- enrichGO( gene = diseasedDGEjoin$entrezgene_id, OrgDb = org.Hs.eg.db, ont = "all", qvalueCutoff = 0.01)

barundifferentiatedGO <- barplot(undifferentiatedGO, showCategory = 15) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=15), 
        legend.key.size = unit(0.6, 'cm'), legend.title = element_text(size=12), legend.text = element_text(size=12)) + 
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave(barundifferentiatedGO, file = paste0(plotDirectory, "undifferentiatedGO.svg"), width = 8, height = 7)

bardifferentiatedGO <- barplot(differentiatedGO, showCategory = 15) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=15), 
        legend.key.size = unit(0.6, 'cm'), legend.title = element_text(size=12), legend.text = element_text(size=12)) + 
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave(bardifferentiatedGO, file = paste0(plotDirectory, 'differentiatedGO.svg'), width = 8, height = 7)

barhealthyGO <- barplot(healthyGO, showCategory = 15) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=15), 
        legend.key.size = unit(0.6, 'cm'), legend.title = element_text(size=12), legend.text = element_text(size=12)) + 
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave(barhealthyGO, file = paste0(plotDirectory, 'healthyGO.svg'), width = 8, height = 7)


bardiseasedGO <- barplot(diseasedGO, showCategory = 15) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=15), 
        legend.key.size = unit(0.6, 'cm'), legend.title = element_text(size=12), legend.text = element_text(size=12)) + 
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave(bardiseasedGO, file = paste0(plotDirectory, 'diseasedGO.svg'), width = 8, height = 7)

#Amount of differentially expressed genes after log2fold filter of abs(1) and padj< 0.01
nrow(topdifferentiatedres.df) #245
nrow(topundifferentiatedres.df) #37
nrow(topdiseasedres.df) #3465
nrow(tophealthyres.df) #3582


topdiffundiffres.df2 <- diffundiffres.df[(diffundiffres.df$baseMean > 10) & (abs(diffundiffres.df$log2FoldChange) > 1), ]
topdiffundiffres.df2 %>% drop_na(log2FoldChange)

topundifferentiatedres.df2 <- undifferentiatedres.df[(undifferentiatedres.df$baseMean > 10) & (abs(undifferentiatedres.df$log2FoldChange) > 1),] 
topundifferentiatedres.df2 %>% drop_na(log2FoldChange)

topdifferentiatedres.df2 <- differentiatedres.df[(differentiatedres.df$baseMean > 10) & (abs(differentiatedres.df$log2FoldChange) > 1), ]
topdifferentiatedres.df2 %>% drop_na(log2FoldChange)

tophealthyres.df2 <- healthyres.df[(healthyres.df$baseMean > 10) & (abs(healthyres.df$log2FoldChange) > 1),]
tophealthyres.df2 %>% drop_na(log2FoldChange)

topdiseasedres.df2 <- diseasedres.df[(diseasedres.df$baseMean) > 10 & (abs(diseasedres.df$log2FoldChange) > 1),]
topdiseasedres.df2 %>% drop_na(log2FoldChange)

nrow(topdiffundiffres.df2) #7278
nrow(topdifferentiatedres.df2) #378
nrow(topundifferentiatedres.df2) #94
nrow(topdiseasedres.df2) #6808
nrow(tophealthyres.df2) #7073









