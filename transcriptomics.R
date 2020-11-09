##steps 1 & 2 of data analysis from a research paper
##1st step: Quality control
library(affy)
raw <- ReadAffy() #read in the files downloaded from the paper
library(affyQCReport)
QCReport(raw, file = "QCReport_for_raw.pdf") #get a Quality control report to help find outliers
library(affyPLM)
fitdata <- fitPLM(raw, normalize = T, background = T) #make data usable for RLE & NUSE
RLE(fitdata) 
NUSE(fitdata) #these two give you boxplots that highlight outliers
library(ggplot2)
matrx <- t(exprs(rma(raw)))
df_pcs <- prcomp(matrx,center = F, scale = F)
df_out <- as.data.frame(df_pcs$x)
group <- factor(c(rep("control",25),rep("cancer",19)))
ggplot(df_out, aes(x=PC1, y = PC2, colour = group, label = as.character(c(1:44)))) + geom_point() + 
  stat_ellipse()+geom_text(aes(label = as.character(c(1:44))), hjust = 0, vjust = 0) #pca plot
##the QCReport, RLE, NUSE, and PCA plot should be compared to identify outlier samples to remove
position <- c(6,8,13,44) 
clean <- raw[-position, ] #these 2 steps are how to remove the samples
rma_clean <- rma(clean, normalize = T, verbose = F)
rma_raw <- rma(raw, normalize = T, verbose = F) #rma() normalizes the datasets
exprs_raw <- exprs(rma_raw)
exprs_clean <- exprs(rma_clean)
boxplot(exprs_raw, col = "brown1") 
boxplot(exprs_clean, col = "aquamarine2") #these boxplots show you why you removed the outliers
##2nd step: Gene annotation
library(hgu133plus2.db)
x <- hgu133plus2.db #package used to exchange PROBEIDs with SYMBOLS
exprs_data <- as.data.frame(exprs_clean)
exprs_data$PROBEID <- as.character(rownames(exprs_data)) 
exprs_annot <- AnnotationDbi::select(x, keys = exprs_data$PROBEID, columns = "SYMBOL")
nonduplicate_anno <- which(duplicated(exprs_annot$PROBEID)==F) #get rid of duplicate PROBEIDs
exprs_annot <- exprs_annot[nonduplicate_anno, ]
exprs_data$SYMBOL <- exprs_annot$SYMBOL
nonduplicate <- which(duplicated(exprs_data$SYMBOL)==F) #get rid of duplicate SYMBOLs too
exprs_data <- exprs_data[nonduplicate, ]
.rowNamesDF(exprs_data, make.names = T) <- exprs_data$SYMBOL #replace PROBEID rownames with SYMBOLS
exprs_data <- subset(exprs_data, select = -c(SYMBOL,PROBEID)) #get rid of the extra columns
##3rd step: Gene filtering & limma
library(limma)
matrix_data <- data.matrix(exprs_data, rownames.force = NA) #make the df a matrix to do quantile()
quantnum <- quantile(matrix_data, 0.04) #assign value to variable for later use
meanofrow <- rowMeans(exprs_data)
exprs_data$meanofrow <- meanofrow
filtered <- exprs_data[exprs_data$meanofrow >=quantnum, ]
filtered <- subset(filtered, select = -meanofrow)
sampletype <- factor(c(rep("control",22),rep("cancer",18))) #factor based off of removed outliers
groupdata <- model.matrix(~factor(sampletype), levels = c("control","cancer"))
fitted_data <- lmFit(filtered,groupdata)
fitted_data <- eBayes(fitted_data)
results_table <- topTable(fitted_data,number = nrow(filtered)) #the final table with all info
##4th step: Making the first set of plots
mm <- as.matrix(filtered)
mmt <- t(filtered)
rownames(mmt)<-sampletype #just to make heatmap easier to read
mmtr <- mmt[ , c(row.names(topt))] #results table except only the top 100 significant genes
library(pheatmap)
pheatmap(mmtr, fontsize_row = 5, fontsize_col = 5,cellwidth = 5, cellheight = 10,
         filename = "heatmap3.pdf") #this heatmap is only of top 100
library(EnhancedVolcano)
EnhancedVolcano(results_table, x = "logFC", y= "adj.P.Val", lab = row.names(results_table),
                col = c('azure3','aquamarine3','cadetblue2','brown2'))
EnhancedVolcano(topt, x = "logFC", y= "adj.P.Val", lab = row.names(topt),
                col = c('azure3','aquamarine3','cadetblue2','brown2'))
##only the first volcano plot is necessary since the 2nd one doesn't show anything interesting
##5th step: Data enrichment and gene ontology plots
results_table$gene_ids <- rownames(results_table)
gvector <- results_table[ , -c(2:6)]
gvectorsorted <- gvector[order(-gvector$logFC), ]
diffexgs <- rownames(gvectorsorted)[abs(gvectorsorted$logFC)>2] #pull out the significant genes
library(clusterProfiler)
library(org.Hs.eg.db)
convert_for_analysis <- bitr(diffexgs, fromType = "SYMBOL", toType = c("ENTREZID","ENSEMBL"),
                             OrgDb = org.Hs.eg.db) 
##ENTREZIDs and ENSEMBLs are needed for the enrichGO() and enrichKEGG() functions
enrichedgo <- enrichGO(gene = convert_for_analysis$ENTREZID, universe = names(diffexgs),
                       OrgDb = org.Hs.eg.db, ont = "CC", readable = T, pAdjustMethod = "fdr")
##the ontology type should put through as "CC", "BP", and "MF" to yield all results
enrichedforplot <- setReadable(enrichedgo, OrgDb = org.Hs.eg.db)
dotplot(enrichedforplot, showCategory = 16) #showCategory won't always show as many as you want
barplot(enrichedforplot, showCategory = 16) #but that is ok, just put the highest possible
plotGOgraph(enrichedforplot) #this plot is hard to read, so you can export it to look at it
eKegg <- enrichKEGG(gene = convert_for_analysis$ENTREZID, organism = "hsa") 
##you need to put "hsa" to let the function know the data is from humans
dotplot(eKegg, showCategory = 3) #the same thing is true for these showCategorys
barplot(eKegg, showCategory = 16)
library(pathview)
pathview(gene.data = eKegg, pathway.id = "hsa04976", species = "hsa")
##pathview should be used on the most significant KEGG pathway, found in the dotplot / barplot


