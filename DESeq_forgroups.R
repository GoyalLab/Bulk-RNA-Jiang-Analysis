#library prep####
library(rhdf5) # For hdf5 files from bootstraps.
library(tidyverse) # Provides helpful R packages for data science.
library(tximport) # How Kallisto results are brought into R.
library(biomaRt) #older versions may use ensembldb instead, below
#library(ensembldb) # Helps interface with ensembl.
library(EnsDb.Hsapiens.v86) # Human-specific database package.
library(beepr) # Surprise function, very necessary.
library(datapasta) #
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
rm(list=ls())
#read in condition data####
healthy <- read_tsv("Healthyonly.tsv")
diseased <- read_tsv("Diseasedonly.tsv")
differentiated <- read_tsv("Differentiatedonly.tsv")
undifferentiated <- read_tsv("Undifferentiatedonly.tsv")

#bring in file paths for abundance data and confirm exist####
differentiatedpath <- file.path(differentiated$...1, "abundance.tsv")
healthypath <- file.path(healthy$...1, "abundance.tsv")
undifferentiatedpath <- file.path(undifferentiated$...1, "abundance.tsv")
diseasedpath <- file.path(diseased$...1, "abundance.tsv")
all(file.exists(differentiatedpath))
all(file.exists(undifferentiatedpath))
all(file.exists(diseasedpath))
all(file.exists(healthypath))

#create Ensembl map in correct form to label data####
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")

#bring in data using file paths and label with gene####
Txi_differentiated <- tximport(differentiatedpath, type = "kallisto", tx2gene = Tx, txOut = FALSE, countsFromAbundance = "no", ignoreTxVersion = TRUE)
Txi_undifferentiated <- tximport(undifferentiatedpath, type = "kallisto", tx2gene = Tx, txOut = FALSE, countsFromAbundance = "no", ignoreTxVersion = TRUE)
Txi_diseased <- tximport(diseasedpath, type = "kallisto", tx2gene = Tx, txOut = FALSE, countsFromAbundance = "no", ignoreTxVersion = TRUE)
Txi_healthy <- tximport(healthypath, type = "kallisto", tx2gene = Tx, txOut = FALSE, countsFromAbundance = "no", ignoreTxVersion = TRUE)

#round count data to prevent decimals so it can be used for analysis####
Txi_differentiated$counts= round(Txi_differentiated$counts)
Txi_undifferentiated$counts= round(Txi_undifferentiated$counts)
Txi_diseased$counts= round(Txi_diseased$counts)
Txi_healthy$counts= round(Txi_healthy$counts)

#separate out just the count data####
differentiatedcounts <- Txi_differentiated$counts
undifferentiatedcounts <- Txi_undifferentiated$counts
diseasedcounts <- Txi_diseased$counts
healthycounts <- Txi_healthy$counts

#label columns for count data####
differentiatedLabels <- differentiated$...1
undifferentiatedLabels <- undifferentiated$...1
diseasedLabels <- diseased$...1
healthyLabels <- healthy$...1

colnames(undifferentiatedcounts) <- c(undifferentiatedLabels)
colnames(diseasedcounts) <- c(diseasedLabels)
colnames(healthycounts) <- c(healthyLabels)
colnames(differentiatedcounts) <- c(differentiatedLabels)

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
all(rownames(differentiated)  == colnames(differentiatedcounts))
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


#diseasedkeep <- rowSums(counts(diseasedddds)) >= 10
#healthykeep <- rowSums(counts(healthyddds)) >= 10
#differentiatedkeep <- rowSums(counts(differentiateddds)) >= 10
#undifferentiatedkeep <- rowSums(counts(undifferentiateddds)) >= 10
#diseasedddds <- diseasedddds[keep,]
#healthyddds <- healthyddds[keep,]
#differentiateddds <- differentiateddds[keep,]
#undifferentiateddds <- undifferentiateddds[keep,]

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
ggplot(plotdifferentiated) + aes(PC1, PC2, colour=Healthy_or_Disease) + geom_text_repel(label = rownames(plotdifferentiated)) + geom_point(size=3) +theme_bw() + labs(title= "Healthy Differentiated Vs Diseased Differentiated", subtitle = "Plot1") + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed() 

#plot2
percentVar <- round(100 * attr(plothealthy, "percentVar"))
ggplot(plothealthy) + aes(PC1, PC2, colour=Differentiated_or_Undifferentiated) + geom_text_repel(label = rownames(plothealthy)) + geom_point(size=3) +theme_bw() + labs(title= "Healthy Differentiated Vs Healthy Undifferentiated", subtitle = "Plot2") + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +  coord_fixed()

#plot3
percentVar <- round(100 * attr(plotundifferentiated, "percentVar"))
ggplot(plotundifferentiated, aes(PC1, PC2, color=Healthy_or_Disease)) + geom_point(size=3) + theme_bw() + geom_text_repel(label = rownames(plotundifferentiated)) + labs(title= "Healthy Undifferentiated Vs Diseased Undifferentiated ", subtitle = "Plot3") + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +  coord_fixed() 

#plot4
percentVar <- round(100 * attr(plotdiseased, "percentVar"))
ggplot(plotdiseased, aes(PC1, PC2, color=Differentiated_or_Undifferentiated)) + geom_point(size=3) + theme_bw() + geom_text_repel(label = rownames(plotdiseased)) + labs(title= "Diseased Differentiated Vs Diseased Undifferentiated", subtitle = "Plot4") + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +  coord_fixed()

#make subset of data####
healthyresOrdered <- healthyres[order(healthyres$pvalue),]
diseasedresOrdered <- diseasedres[order(diseasedres$pvalue),]
differentiatedresOrdered <- differentiatedres[order(differentiatedres$pvalue),]
undifferentiatedresOrdered <- undifferentiatedres[order(undifferentiatedres$pvalue),]

#first filter
healthyresSig <- subset(healthyresOrdered, padj < 0.05)
diseasedresSig <- subset(diseasedresOrdered, padj < 0.05)
differentiatedresSig <- subset(differentiatedresOrdered, padj < 0.05)
undifferentiatedresSig <- subset(undifferentiatedresOrdered, padj < 0.05)

#write.csv(as.data.frame(healthyresSig), 
#          file="Healthydifferentiated_Healthyundifferentiated.csv")

#write.csv(as.data.frame(diseasedresSig), 
#          file="Diseaseddifferentiated_Diseasedundifferentiated.csv")

#write.csv(as.data.frame(differentiatedresSig), 
#          file="healthyDifferentiated_diseasedDifferentiated.csv")

#write.csv(as.data.frame(undifferentiatedresSig), 
#          file="healthyUnifferentiated_diseasedUnifferentiated.csv")


#healthylog2_changes <- healthyresSig[, "log2FoldChange"]
#diseasedlog2_changes <- diseasedresSig[, "log2FoldChange"]
#differentiatedlog2_changes <- differentiatedresSig[, "log2FoldChange"]
#undifferentiatedlog2_changes <- undifferentiatedresSig[, "log2FoldChange"]

#turn into data frame
undifferentiatedres.df <- as.data.frame(undifferentiatedresSig)
differentiatedres.df <- as.data.frame(differentiatedresSig)
healthyres.df <- as.data.frame(healthyresSig)
diseasedres.df <- as.data.frame(diseasedresSig)

#second filter, removing insignificant data and genes with NA values for log2 fold change
topundifferentiatedres.df <- undifferentiatedres.df[(undifferentiatedres.df$baseMean > 10) & (abs(undifferentiatedres.df$log2FoldChange) > 1),] 
topundifferentiatedres.df %>% drop_na(log2FoldChange)

topdifferentiatedres.df <- differentiatedres.df[(differentiatedres.df$baseMean > 10) & (abs(differentiatedres.df$log2FoldChange) > 1), ]
topdifferentiatedres.df %>% drop_na(log2FoldChange)
                             
tophealthyres.df <- healthyres.df[(healthyres.df$baseMean > 10) & (abs(healthyres.df$log2FoldChange) > 1),]
tophealthyres.df %>% drop_na(log2FoldChange)

topdiseasedres.df <- diseasedres.df[(diseasedres.df$baseMean) > 10 & (abs(diseasedres.df$log2FoldChange) > 1),]
topdiseasedres.df %>% drop_na(log2FoldChange)

#ordering the data to have the top log2fold genes on top and the lowest non NA values at the bottom
topundifferentiatedres.df <- topundifferentiatedres.df[order(topundifferentiatedres.df$log2FoldChange, decreasing = TRUE),]
topdifferentiatedres.df <- topdifferentiatedres.df[order(topdifferentiatedres.df$log2FoldChange, decreasing = TRUE),]
tophealthyres.df <- tophealthyres.df[order(tophealthyres.df$log2FoldChange, decreasing = TRUE),]
topdiseasedres.df <- topdiseasedres.df[order(topdiseasedres.df$log2FoldChange, decreasing = TRUE),]

#write.csv(as.data.frame(tophealthyres.df), 
#          file="FilteredHealthydifferentiated_Healthyundifferentiated.csv")

#write.csv(as.data.frame(topdiseasedres.df), 
#          file="FilteredDiseaseddifferentiated_Diseasedundifferentiated.csv")

#write.csv(as.data.frame(topdifferentiatedres.df), 
#          file="FilteredhealthyDifferentiated_diseasedDifferentiated.csv")

#write.csv(as.data.frame(topundifferentiatedres.df), 
#          file="FilteredhealthyUnifferentiated_diseasedUnifferentiated.csv")

#create a gene row so its easier to search
topundifferentiatedres.df$gene <- rownames(topundifferentiatedres.df)
topdifferentiatedres.df$gene <- rownames(topdifferentiatedres.df)
topdiseasedres.df$gene <- rownames(topdiseasedres.df)
tophealthyres.df$gene <- rownames(tophealthyres.df)

# since the data is organised by log numbers we can find the row numbers of the 50 most up and down regulated genes
num_keep <- 50
topdifferentiatedrowkeep <- c(seq(1:num_keep), seq((nrow(topdifferentiatedres.df)-num_keep), nrow(topdifferentiatedres.df)))
topundifferentiatedrowkeep <- c(seq(1:num_keep), seq((nrow(topundifferentiatedres.df)-num_keep), nrow(topundifferentiatedres.df)))
topdiseasedrowkeep <- c(seq(1:num_keep), seq((nrow(topdiseasedres.df)-num_keep), nrow(topdiseasedres.df)))
tophealthyrowkeep <- c(seq(1:num_keep), seq((nrow(tophealthyres.df)-num_keep), nrow(tophealthyres.df)))


#filteredhealthylog2_changes <- as.matrix(tophealthyres.df[tophealthyrowkeep,]$log2FoldChange)
#colnames(filteredhealthylog2_changes) <- "Healthy Differentiated Vs Healthy Undifferentiated Log2- fold change"

#filtereddiseasedlog2_changes <- as.matrix(topdiseasedres.df[topdiseasedrowkeep,]$log2FoldChange)
#colnames(filtereddiseasedlog2_changes) <- "Diseased Differentiated Vs Diseased Undifferentiated Log2- fold change"

#filtereddifferentiatedlog2_changes <- as.matrix(topdifferentiatedres.df[topdifferentiatedrowkeep,]$log2FoldChange)
#colnames(filtereddifferentiatedlog2_changes) <- "Healthy Differentiated Vs Diseased Differentiated Log2- fold change"

#filteredundifferentiatedlog2_changes <- as.matrix(topundifferentiatedres.df[topundifferentiatedrowkeep,]$log2FoldChange)
#colnames(filteredundifferentiatedlog2_changes) <- "Healthy Undifferentiated Vs Diseased Undifferentiated Log2- fold change"

topundifferentiatedres.df$gene <- rownames(topundifferentiatedres.df)
topdifferentiatedres.df$gene <- rownames(topdifferentiatedres.df)
topdiseasedres.df$gene <- rownames(topdiseasedres.df)
tophealthyres.df$gene <- rownames(tophealthyres.df)
 
# merging the count data will filter out some insignificant  genes
allundifferentiated <- merge(undifferentiatedcounts, topundifferentiatedres.df, by = 0)
alldifferentiated <- merge(differentiatedcounts, topdifferentiatedres.df, by = 0)
alldiseased <- merge(diseasedcounts, topdiseasedres.df, by = 0)
allhealthy <- merge(healthycounts, tophealthyres.df, by = 0)

#the log data was useful for filtering but the only information to make the heat map is count data, sample name and gene name. Remove anything else.
allundifferentiatedcounts <- allundifferentiated[,2:11]
alldifferentiatedcounts <- alldifferentiated[,2:11]
allhealthycounts <- allhealthy[,2:11]
alldiseasedcounts <- alldiseased[,2:11]

#label genes
rownames(allundifferentiatedcounts) <- topundifferentiatedres.df$gene
rownames(alldifferentiatedcounts) <- topdifferentiatedres.df$gene
rownames(alldiseasedcounts) <- topdiseasedres.df$gene
rownames(allhealthycounts) <- tophealthyres.df$gene

#create a data frame that only has the most up or down regulated genes
filteredallhealthycounts <- as.data.frame(allhealthycounts[tophealthyrowkeep,])
filteredalldiseasedcounts <- as.data.frame(alldiseasedcounts[topdiseasedrowkeep,])
filteredalldifferentiatedcounts <- as.data.frame(alldifferentiatedcounts[topdifferentiatedrowkeep,])
filteredallundifferentiatedcounts <- as.data.frame(allundifferentiatedcounts[topundifferentiatedrowkeep,])

#make heat maps using our filtered out 100 genes (50 most unregulated, 50 most down regulated)####
library(pheatmap)
#heatmap1
pheatmap(filteredallundifferentiatedcounts, scale = "row", main = "Undifferentiated Samples", fontsize_row = 7)

#heatmap2
pheatmap(filteredalldifferentiatedcounts, scale = "row", main = "Differentiated Samples", fontsize_row = 7)

#heatmap3
pheatmap(filteredalldiseasedcounts, scale = "row", main = "Diseased Samples", fontsize_row = 7)

#heatmap4
pheatmap(filteredallhealthycounts, scale = "row", main = "Healthy Samples", fontsize_row = 7)

#filter data again to make volcano plots####
library(ggplot2)
library(ggrepel)

mycolours <- c("blue", "red", "black")
names(mycolours) <- c("Down regulated", "Up regulated", "Nonsignificant")

#order of pvalue doesnt matter
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
subhealthyresOrdered$differentialexpression[subhealthyresOrdered$log2FoldChange > 1 & subhealthyresOrdered$padj < 0.05] <- "Up regulated"
subhealthyresOrdered$differentialexpression[subhealthyresOrdered$log2FoldChange < -1 & subhealthyresOrdered$padj < 0.05] <- "Down regulated"
subhealthyresOrdered$deLabel <- NA
subhealthyresOrdered$deLabel[subhealthyresOrdered$differentialexpression != "Nonsignificant"] <- rownames(subhealthyresOrdered)[subhealthyresOrdered$differentialexpression != "Nonsignificant"]
#volcano1
ggplot(subhealthyresOrdered, aes(x= log2FoldChange, y= -log10(padj), col = differentialexpression, label = deLabel)) + geom_point() + geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed") + geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") + geom_text_repel() + scale_color_manual(values = mycolours) + labs(title= "Healthy Differentiated Vs Healthy Undifferentiated") + xlab(bquote(log[2]*FoldChange)) + ylab(bquote(-log[10]*(padj)))


subdiseasedresOrdered$differentialexpression <- "Nonsignificant"
subdiseasedresOrdered$differentialexpression[subdiseasedresOrdered$log2FoldChange > 1 & subdiseasedresOrdered$padj < 0.05] <- "Up regulated"
subdiseasedresOrdered$differentialexpression[subdiseasedresOrdered$log2FoldChange < -1 & subdiseasedresOrdered$padj < 0.05] <- "Down regulated"
subdiseasedresOrdered$deLabel <- NA
subdiseasedresOrdered$deLabel[subdiseasedresOrdered$differentialexpression != "Nonsignificant"] <- rownames(subdiseasedresOrdered)[subdiseasedresOrdered$differentialexpression != "Nonsignificant"]
#volcano2
ggplot(subdiseasedresOrdered, aes(x= log2FoldChange, y= -log10(padj), col = differentialexpression, label = deLabel)) + geom_point() + geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed") + geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") + geom_text_repel() + scale_color_manual(values = mycolours)+ labs(title= "Healthy Undifferentiated Vs Diseased Undifferentiated ") + xlab(bquote(log[2]*FoldChange)) + ylab(bquote(-log[10]*(padj)))


subdifferentiatedresOrdered$differentialexpression <- "Nonsignificant"
subdifferentiatedresOrdered$differentialexpression[subdifferentiatedresOrdered$log2FoldChange > 1 & subdifferentiatedresOrdered$padj < 0.05] <- "Up regulated"
subdifferentiatedresOrdered$differentialexpression[subdifferentiatedresOrdered$log2FoldChange < -1 & subdifferentiatedresOrdered$padj < 0.05] <- "Down regulated"
subdifferentiatedresOrdered$deLabel <- NA
subdifferentiatedresOrdered$deLabel[subdifferentiatedresOrdered$differentialexpression != "Nonsignificant"] <- rownames(subdifferentiatedresOrdered)[subdifferentiatedresOrdered$differentialexpression != "Nonsignificant"]
#volcano3
ggplot(subdifferentiatedresOrdered, aes(x= log2FoldChange, y= -log10(padj), col = differentialexpression, label = deLabel)) + geom_point() + geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed") + geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") + geom_text_repel() + scale_color_manual(values = mycolours) + labs(title= "Diseased Differentiated Vs Diseased Undifferentiated") + xlab(bquote(log[2]*FoldChange)) + ylab(bquote(-log[10]*(padj)))


subundifferentiatedresOrdered$differentialexpression <- "Nonsignificant"
subundifferentiatedresOrdered$differentialexpression[subundifferentiatedresOrdered$log2FoldChange > 1 & subundifferentiatedresOrdered$padj < 0.05] <- "Up regulated"
subundifferentiatedresOrdered$differentialexpression[subundifferentiatedresOrdered$log2FoldChange < -1 & subundifferentiatedresOrdered$padj < 0.05] <- "Down regulated"
subundifferentiatedresOrdered$deLabel <- NA
subundifferentiatedresOrdered$deLabel[subundifferentiatedresOrdered$differentialexpression != "Nonsignificant"] <- rownames(subundifferentiatedresOrdered)[subundifferentiatedresOrdered$differentialexpression != "Nonsignificant"]
#volcano4
ggplot(subundifferentiatedresOrdered, aes(x= log2FoldChange, y= -log10(padj), col = differentialexpression, label = deLabel)) + geom_point() + geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed") + geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") + geom_text_repel() + scale_color_manual(values = mycolours) + labs(title= "Healthy Undifferentiated Vs Diseased Undifferentiated ") + xlab(bquote(log[2]*FoldChange)) + ylab(bquote(-log[10]*(padj)))






